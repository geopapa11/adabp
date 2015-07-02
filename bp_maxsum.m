function [map_seq, map_val, b, m, m_idx, nghb] = bp_maxsum(E, phi, psi, obs, varargin)
% function [map_seq, map_val, b, m, m_idx, nghb] = bp_maxsum(E, phi, psi, obs, varargin)
% 
% This function performs max-sum algorithm in a forest, in other words, in
% a graph whose connected components are trees.
% It first propagates all the messages from the leaves all the way up to
% the root. Then, it retrieves the MAP sequence by backtracking from the
% root down to the leaves.
%
% INPUT
%        N:  1x1   scalar        -  Number of nodes (variables)
%       Ne:  1x1   scalar        -  Number of edges
%     card:  Nx1   double array  -  Cardinality (alphabet) for each variable
%                                   That is, each variable x_i takes values
%                                   in {1, ..., card(i)}
%        E:  Nex2 double matrix  -  Edges
%                                   Each element corresponds to an edge
%                                   (pair of nodes that are linked).
%      phi:  Nx1 cell array      -  Node potentials
%            or                     Each element is a card(i)x1 vector
%            card x 1               (if phi is common among all nodes)
%            double array
%      psi:  Nex1 cell array     -  Pairwise potentials
%                                   Each element represents a pairwise 
%                                   potential between nodes defined by the 
%                                   corresponding edge in matrix E.
%            or                     Each element is a card(C(i,1))xcard(C(i,2)) vector
%            card x card            (if psi is common among all edges)
%            double matrix
%      obs:  1xN double array    -  Observed values
%                                   Each element equals to the value of the 
%                                   variable in the corresponding order.
%                                   For instance, if three variables have
%                                   labels 1, 2, 3, the 1st, 2nd and 3rd
%                                   elements correspond to variables 1, 2 and 3.
%                                   If three variables have labels 4, 7, 9,
%                                   the 1st, 2nd and 3rd elements correspond
%                                   to variables 4, 7 and 9.
%                                   (If a variable is unobserved, the
%                                   corresponding element is set to NaN)
%  varargin: cell array          -  Contains pairs of arguments
%                                   * If varargin{i} = 'root', then varargin{i+1}
%                                     is an 1x1 scalar representing the root.
%                                     (Default value: root chosen randomly)
% 
% 
% OUTPUT
%  map_seq:  1xN double matrix  -  MAP sequence
%                                  Sequence of variable values corresponding
%                                  to the maximum value of this function
%  map_val:  1x1 scalar         -  MAP value
%                                  Maximum value of this function
%        b:  Nx1 cell array     -  Node max-marginal
%                                  Each element is: b{i} = max_{x_j | j ~= i} p(x_1,x_2,...,x_N)
%        m:  Nx1 cell array     -  Messages
%                                  Each row i is a cell array representing 
%                                  the messages that are incoming to node i
%    m_idx:  Nx1 cell array     -  Message maximum indices
%                                  Each row lists the indices of the source
%                                  nodes that maximize the value of the
%                                  message at the target (i) node
%     nghb:  Nx1 or 1xN         -  Neighbors
%            cell array            Each element lists the neighbors of the
%                                  corresponding variable
% 
% 
% EXAMPLES
%
% We show an example of max-product in a graph whose connected components
% are trees.
% 
%       1             
%     /   \
%    2     3        5 -- 6
%           \
%            4
% 
% Assume we have the following edges, node and pairwise potentials:
% 
% E = [1  2;
%      1  3;
%      5  6;
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.4 .6]';  phi{3} = [.9 .1]';  phi{4} = [.1 .9]';
% phi{5} = [.4 .6]';  phi{6} = [.9 .1]';
% 
% psi{1} = [.9 .1; 
%           .3 .7];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.8 .2; 
%           .1 .9];
% psi{4} = [.2 .8;
%           .7 .3];
% 
% [map_seq, map_val, b, m, m_idx, nghb] = bp_maxsum(E, phi, psi);
% 
% This will produce output:
% map_seq = [];
% map_val = ;
% b{1} = [0.5470  0.4530]';  b{2} = [0.5695  0.4305]';
% b{3} = [0.8923  0.1077]';  b{4} = [0.0463  0.9537]';
% b{5} = [0.7327  0.2673]';  b{6} = [0.8465  0.1535]';
% 
% Author: geopapa
% $ Date: 2014/05/15 15:25:12 $

    N       = length(phi);
    Ne      = size(E,1);
    nghb    = cell(N,1);
    map_seq = zeros(1,N);
    m       = cell(N,1);
    m_idx   = cell(N,1);
    b       = cell(N,1);
    
    % Make sure all node potentials are in vector form
    Icol       = cell2mat(cellfun(@iscolumn, phi, 'UniformOutput', false));
    phi(~Icol) = cellfun(@transpose,  phi(~Icol), 'UniformOutput', false);
    
    % Determine the root
    if ~any(strcmp(varargin,'root'))
        root     = randi(N);  % Choose a root randomly
        root_idx = [];
    else
        root     = varargin{find(strcmp(varargin,'root'))+1};
        root_idx = find(E==root,1);
        if isempty(root_idx),  error('You entered a root that does not participate in any edge.');  end
    end
    
    % Relabel nodes from 1 to N
    if ~isempty(E),         E    = relabel(E);   end
    if ~isempty(root_idx),  root = E(root_idx);  end
    
    % Find the neighbors of each node
    for i = 1:N
        [row,col]   = find(E==i);
        col         = col+1;
        col(col==3) = 1;
        idx         = sub2ind(size(E),row,col);
        nghb{i}     = E(idx)';
    end
    nghb = cellfun(@sort, nghb, 'UniformOutput', false);
    
    % Get the connected components of the graph
    comp = conn_comp(nghb);
    
    % Each element of edge2comp gives the component that an edge belongs to
    edge2comp = zeros(Ne,1);
    if ~isempty(E)
        for i = 1:length(comp)
            I = ismember(E(:,1), comp{i});
            edge2comp(I) = i;
        end
    end
    
    % Run bp for each connected component
    if length(comp)==1
        if nargin >= 4 && ~isempty(obs)
            [map_seq, map_val, b, m, m_idx] = bp_maxsum_core(E, phi, psi, obs, 'root', root, 'nghb', nghb);
        else
            [map_seq, map_val, b, m, m_idx] = bp_maxsum_core(E, phi, psi,  [], 'root', root, 'nghb', nghb);
        end
    else
        map_val = 1;
        for i = 1:length(comp)
            E_cc    = E(edge2comp==i,:);
            phi_cc  = phi(comp{i});
            psi_cc  = psi(edge2comp==i);
            nghb_cc = nghb(comp{i});
            if nargin >= 4 && ~isempty(obs)
                obs_cc       = obs(comp{i});
                [map_seq_cc, map_val_cc, b_cc, m_cc, m_idx_cc] = bp_maxsum_core(E_cc, phi_cc, psi_cc, obs_cc, 'nghb', nghb_cc);
            else
                [map_seq_cc, map_val_cc, b_cc, m_cc, m_idx_cc] = bp_maxsum_core(E_cc, phi_cc, psi_cc, 'nghb', nghb_cc);
            end
            map_seq(comp{i}) = map_seq_cc;
            map_val          = map_val*map_val_cc;
            b(comp{i})       = b_cc;
            m(comp{i})       = m_cc;
            m_idx(comp{i})   = m_idx_cc;
        end
    end
end