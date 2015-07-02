function [map_seq, map_val, b, msg, msg_idx, nghb] = bp_maxsum_core(E, phi, psi, obs, varargin)
% function [map_seq, map_val, b, msg, msg_idx, nghb] = bp_maxsum_core(E, phi, psi, obs, varargin)
% 
% This function performs max-sum algorithm in a tree graph, where only 
% pairwise cliques exist.
% It first propagates all the messages from the leaves all the way up to
% the root. Then, it retrieves the MAP sequence by backtracking from the
% root down to the leaves.
% 
% Algorithm:
% 
% 1. Initialization
% 
%    for i=1:N
%       for j=1:N(i)
%           m_{i->j}^{0}(x_j) = 0
%       end
%    end
% 
% 
% 2. Apply from children to parents (start from leaves, end at the root)
%    blf_{i}(x_i)   = log(phi(x_i)) + sum_{k \in N(i)} m_{k->i}(x_i)
%    m_{i->j}(x_j)  = max_{x_i}    { log(psi(x_i,x_j) + blf_{i}(x_i) - m_{j->i}(x_i) }
%   idx_{i->j}(x_j) = argmax_{x_i} { log(psi(x_i,x_j) + blf_{i}(x_i) - m_{j->i}(x_i) },
%   where j is the parent of i: j = pa(i).
% 
% 
% 3. Do the same from root down to leaves
% 
% 
% 4. Compute the root's (x_r) max-marginal
%       b_r(x_r) = log(phi(x_r)) + sum_{k \in N(r)} m_{k->r}(x_r)
% 
% 
% 5. Backtrack to recover a MAP configuration
%    
%    Start from the root
%       x_r^{map} = argmax_{x_r} { b_r(x_r) }
% 
%    Continue down at the leaves
%       for j \in N(i)
%           x_j^{map} = idx_{j->i}(x_i^{map})
%       end
%   Continue backtracking until all nodes' maximum values have been tracked.
% 
% 
% 6. Log-Likelihood and likelihood of maximum configuration are
% 
%    log(p(x_1^{map}, ..., x_N^{map})) = b_i{x_i^{map}
%    p(x_1^{map}, ..., x_N^{map})      = exp(b_i{x_i^{map})  ,
%    respectively.
%
% 
% INPUT
%          N:   1x1   scalar         -  Number of nodes (variables)
%         Ne:   1x1   scalar         -  Number of edges
%       card:   1x1   scalar         -  Cardinality (alphabet) of the variables
%                                       That is, each variable x_i takes values 
%                                       from {1, ..., card}, that is, 
%                                       X_i \in \mathcal{X}_i = {1, ..., card_i}
%          E: Nex2 double matrix     -  Pairwise cliques
%                                       Each row corresponds to an edge
%        phi:  Nx1 cell array        -  Node potentials
%                                       (each element is a |X_i| x 1 vector)
%         or  card x 1 double array  -  Node potential 
%                                       (common across all nodes)
%        psi:  Nex1 cell array       -  Pairwise potentials
%                                       Each element represents a pairwise potential 
%                                       between nodes defined by double matrix E
%                                       (each element (i,j) is a |X_i|x|X_j| vector)
%         or  card x card            -  Pairwise potential 
%             double matrix             (common across all edges)
%       obs:  1xN double array       -  Observed values
%                                       Each element equals to the value of the 
%                                       variable in the corresponding order.
%                                       For instance, if three variables have
%                                       labels 1, 2, 3, the 1st, 2nd and 3rd
%                                       elements correspond to variables 1, 2 and 3.
%                                       If three variables have labels 4, 7, 9,
%                                       the 1st, 2nd and 3rd elements correspond
%                                       to variables 4, 7 and 9.
%                                       (If a variable is unobserved, the
%                                       corresponding element is set to NaN
%                                       In addition, if obs=[], it means no
%                                       values are observed)
%  varargin: cell array              -  Contains pairs of arguments
%                                       * If varargin{i} = 'root', then varargin{i+1}
%                                         is an 1x1 scalar representing the root.
%                                         (Default value: root chosen randomly)
%                                       * If varargin{i} = 'nghb', then varargin{i+1}
%                                         are the neighbors structure of
%                                         this graph
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
%      msg:  Nx1 cell array     -  Messages
%                                  Each row i is a cell array representing 
%                                  the messages that are incoming to node i
%  msg_idx:  Nx1 cell array     -  Message maximum indices
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
% E.g. assume the simple graphical model:
%
%       1
%     /   \
%    2     3
%           \
%            4
% 
% 
% We have the following edges, node and pairwise potentials, respectively:
% 
% E = [1  2;
%      1  3; 
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.2 .3 .4]';  phi{3} = [.9 .1]';  phi{4} = [.2 .1 .6 .1]';
% 
% psi{1} = [.1 .7 .2; 
%           .3 .6 .1];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.1 .2 .65 .05;
%           .1 .7 .05 .15];
%
% [map_seq, map_val, b, msg, msg_idx, nghb] = bp_maxsum_core(E, phi, psi);
% 
% This will produce output:
% map_seq = [1  2  1  3];
% 
% 
% Lastly, we will make an example with some variables being observed.
% 
% E = [1  2;
%      1  3; 
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.2 .3 .4]';  phi{3} = [.9 .1]';  phi{4} = [.2 .1 .6 .1]';
% 
% psi{1} = [.1 .7 .2; 
%           .3 .6 .1];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.1 .2 .65 .05;
%           .1 .7 .05 .15];
% 
% obs    = [NaN 1 NaN 2]';  % Remember NaN indicates an unobserved variable
% 
% [map_seq, map_val, b, msg, msg_idx, nghb] = bp_maxsum_core(E, phi, psi, obs);
% 
% This will produce output:
% map_seq = [2  1  1  2];
%
% Author: geopapa
% $ Date: 2014/05/15 00:21:33 $

    % Check if there is only one node in the graph
    if isempty(E)
        if ~iscell(phi)
            [map_val, map_seq] = max(phi);
        elseif iscell(phi) && length(phi)==1
            [map_val, map_seq] = max(phi{1});
        else
            error('Edge matrix is empty and there are more than one nodes in the graph');
        end
        b       = [];
        msg     = [];
        msg_idx = [];
        nghb    = [];
        return;
    end

    if size(E,1)~=length(psi) && iscell(psi)
        error('The number of edges should equal the number of pairwise potentials.');
    end
    
    % If node potential is common among all nodes
    if ~iscell(phi)
        tmp     = phi;
        unq_nds = unique(E);
        phi     = cell(1,length(unq_nds));
        phi(:)  = {tmp};
        clear tmp;
    end

    N    = length(phi);           % number of nodes
    card = cellfun(@length,phi);  % cardinality (alphabet) for each node
    b    = cell(N,1);
        
    % Make sure all node potentials are in vector form    
    Icol       = cell2mat(cellfun(@iscolumn, phi, 'UniformOutput', false));
    phi(~Icol) = cellfun(@transpose, phi(~Icol),  'UniformOutput', false);
    
    % Check if there are observed variables
    if nargin>=4 && ~isempty(obs)
        if ~isrow(obs)   % convert to row vector
            obs = obs';
        end
        
        if any(obs<1) || any(obs>card)
            error('There are observed values that are outside their specified alphabet.');
        end
        
        idx = find(~isnan(obs));
        for i = 1:length(idx)
            I = true(card(idx(i)),1);
            I(obs(idx(i))) = 0;
            phi{idx(i)}(I) = 0;
        end
    end
    
    if ~isempty(E),  E = relabel(E);  end   % Relabel nodes from 1 to N
        
    % Find the neighbors of each node
    if any(strcmp(varargin,'nghb'))
        nghb = varargin{find(strcmp(varargin,'nghb'))+1};
        nghb = relabel(nghb);
    else
        nghb = cell(N,1);
        for i = 1:N
            [row,col]   = find(E==i);
            col         = col+1;
            col(col==3) = 1;
            idx         = sub2ind(size(E),row,col);
            nghb{i}     = E(idx)';
        end
        nghb = cellfun(@sort, nghb, 'UniformOutput', false);
    end
    
    % Create a map of edges to their linear indices in matrix E
    keys     = unique(E(:,1));
    edgeIdx  = (1:size(E,1))';
    edgeOrig = sparse(keys,1,1:length(keys),max(E(:)),1);
    edge2idx = cell(length(keys),1);
    for k = 1:length(edge2idx)
        I           = E(:,1)==keys(k);
        edge2idx{k} = sparse(E(I,2),1,edgeIdx(I),max(E(I,2)),1);
    end
    
    % Create a map of neighbors to their linear indices in cell matrix nghb
    nghb2idx = cell(N,1);
    for i = 1:N
        keys        = nghb{i};
        vals        = 1:length(nghb{i});
        nghb2idx{i} = sparse(keys,1,vals,max(keys),1);
    end
    
    % Determine the root
    if ~any(strcmp(varargin,'root'))
        root = randi(N);  % Choose a root randomly
    else
        root = varargin{find(strcmp(varargin,'root'))+1};
    end
    
    % Find the children of each node
    ch = get_ch(root, nghb);
    
    % Generate the message schedule
    queue = msg_schedule(root, ch);
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Initialize the messages  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg = cell(N,1);
    for i = 1:N
        tmp    = cell(size(nghb{i}));
        tmp(:) = {zeros(card(i),1)};
        msg{i} = tmp;
    end
    msg_idx  = msg;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  1st Phase: Leaves -> Root  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for idx = size(queue,1):-1:1
        orig = queue(idx,1);
        targ = queue(idx,2);
        
        blf = log(phi{orig}) + sum(cell2mat(msg{orig}),2);  % Consider all the incoming messages to node orig
        
        if iscell(psi)
            if edgeOrig(orig) ~= 0 && edge2idx{edgeOrig(orig)}(targ) ~= 0
                psi_cur = psi{edge2idx{edgeOrig(orig)}(targ)};
            else
                psi_cur = psi{edge2idx{edgeOrig(targ)}(orig)}';
            end
        else
            if edgeOrig(orig) ~= 0 && edge2idx{edgeOrig(orig)}(targ) ~= 0
                psi_cur = psi;
            else
                psi_cur = psi';
            end
        end
                
        [max_val, max_idx]                  = max(bsxfun(@plus, log(psi_cur), blf - msg{orig}{nghb2idx{orig}(targ)}));  % subtract m_{targ->orig}
        msg{targ}{nghb2idx{targ}(orig)}     = max_val(:);  % m_{targ->orig}
        msg_idx{targ}{nghb2idx{targ}(orig)} = max_idx(:);  % indices of that maximize x_{orig} for every value of x_{targ}
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  2nd Phase: Root -> Leaves  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for idx = 1:size(queue,1)
        orig = queue(idx,2);
        targ = queue(idx,1);
        
        blf = log(phi{orig}) + sum(cell2mat(msg{orig}),2);  % Consider all the incoming messages to node orig
        
        if iscell(psi)
            if edgeOrig(orig) ~= 0 && edge2idx{edgeOrig(orig)}(targ) ~= 0
                psi_cur = psi{edge2idx{edgeOrig(orig)}(targ)};
            else
                psi_cur = psi{edge2idx{edgeOrig(targ)}(orig)}';
            end
        else
            if edgeOrig(orig) ~= 0 && edge2idx{edgeOrig(orig)}(targ) ~= 0
                psi_cur = psi;
            else
                psi_cur = psi';
            end
        end
        
        [max_val, max_idx]                  = max(bsxfun(@plus, log(psi_cur), blf - msg{orig}{nghb2idx{orig}(targ)}));  % subtract m_{targ->orig}
        msg{targ}{nghb2idx{targ}(orig)}     = max_val(:);  % m_{targ->orig}
        msg_idx{targ}{nghb2idx{targ}(orig)} = max_idx(:);  % indices of that maximize x_{orig} for every value of x_{targ}
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%  Max-Marginals  %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N        
        b{i} = log(phi{i}) + sum(cell2mat(msg{i}),2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%  Backtracking %%%  (to retrieve the max sequence)
    %%%%%%%%%%%%%%%%%%%%%
    
    % Find x_{root}^{map}
    [~, mroot_idx] = max(b{root});
    
    % Backtracking
    map_seq = btrack(root, mroot_idx, nghb, msg_idx, 'nghb2idx', nghb2idx);
        
    % Evaluate log-likelihood of MAP sequence
    lg_map_val = loglik(map_seq, E, phi, psi);
    map_val    = exp(lg_map_val);
end