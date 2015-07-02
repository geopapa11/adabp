function [b, msg, nghb] = bp_core(E, phi, psi, obs, varargin)
% [b, msg, nghb] = bp_core(E, phi, psi, obs, varargin)
% 
% This function performs belief propagation in a tree graph where only 
% pairwise cliques exist by first propagating all the messages from the
% leaves all the way up to the root and then backwards.
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
%                                   * If varargin{i} = 'nghb', then varargin{i+1}
%                                     are the neighbors structure of this 
%                                     graph
% 
% 
% OUTPUT
%        b:  Nx1 cell array    -  Node marginal
%                                 Each element is: b{i} = p(x_i)
%      msg:  Nx1 cell matrix   -  Messages
%                                 Each row i is a cell array representing 
%                                 the messages that are incoming to node i
%     nghb:  Nx1 or 1xN        -  Neighbors
%            cell array           Each element lists the neighbors of the
%                                 corresponding variable
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
% with edges, node and pairwise potentials:
% 
% E = [1  2;
%      1  3;
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.4 .6]';  phi{3} = [.9 .1]';  phi{4} = [.1 .9]';
% 
% psi{1} = [.9 .1; 
%           .3 .7];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.8 .2; 
%           .1 .9];
% 
% [b, msg] = bp_core(E, phi, psi);
% 
% This will produce output:
% b{1} = [0.6688  0.3312]';  b{2} = [0.6468  0.3532];
% b{3} = [0.5470  0.4530]';  b{4} = [0.1738  0.8262];
% 
% 
% Lastly, we will make an example with some variables being observed.
% The observed nodes are x_2 and x_4 with values 1 and 2, respectively.
% Therefore: obs = [NaN 1 NaN 2]';  % Remember NaN indicates an unobserved
% variable.
% We have the following edges, node and pairwise potentials:
% 
% E = [1  2;
%      1  3;
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.4 .6]';  phi{3} = [.9 .1]';  phi{4} = [.1 .9]';
% 
% psi{1} = [.9 .1; 
%           .3 .7];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.8 .2; 
%           .1 .9];
% 
% obs    = [NaN 1 NaN 2]';
% 
% [b, msg] = bp_core(E, phi, psi, obs);
% 
% This will produce output:
% b{1} = [0.9  0.1]';  b{2} = [1  0];  b{3} = [0.375  0.625]';  b{4} = [0  1]';
% 
% Author: geopapa
% $ Date: 2013/10/15 14:01:36 $

    % Check if there is only one node in the graph
    if isempty(E)
        if ~iscell(phi)
            b = {phi/sum(phi)};
        elseif iscell(phi) && length(phi)==1
            b = {phi{1}/sum(phi{1})};
        else
            error('Edge matrix is empty and there are more than one nodes in the graph');
        end
        msg  = {{ones(size(b{1}))}};
        nghb = [];
        return;
    end
    
    if size(E,1)~=length(psi) && iscell(psi)
        error('The number of edges should equal the number of pairwise potentials.');
    end
    
    % If node potential is common among all nodes, create a separate one for each node
    if ~iscell(phi)
        tmp     = phi;
        unq_nds = unique(E);
        phi     = cell(1,length(unq_nds));
        phi(:)  = {tmp};
        clear tmp;
    end
    
    % Make sure all node potentials are in vector form
    Icol       = cell2mat(cellfun(@iscolumn, phi, 'UniformOutput', false));
    phi(~Icol) = cellfun(@transpose, phi(~Icol),  'UniformOutput', false);
    
    N    = length(phi);           % number of nodes
    Ne   = size(E,1);             % number of edges
    card = cellfun(@length,phi);  % cardinality (alphabet) for each node
    
    % Check if there are observed variables
    if nargin >= 4 && ~isempty(obs)
        % Convert to row vector
        if ~isrow(obs) 
            obs = obs';
        end
        
        if any(obs<1) || any(obs>card)
            error('There are observed values that are outside their specified alphabet.');
        end
        
        idx = find(~isnan(obs));
        for i = 1:length(idx)
            I              = true(card(idx(i)),1);
            I(obs(idx(i))) = 0;
            phi{idx(i)}(I) = 0;
        end        
    end
    
    % Scale up the node and edge potentials if they are too small
    for i = 1:N
        if sum(phi{i}) < 1
            phi{i} = phi{i}/sum(phi{i});
        end                
    end
    
    if iscell(psi)
        for e = 1:Ne
            if sum(psi{e}(:)) < 1
                psi{e} = psi{e}/sum(psi{e}(:));
            end
        end
    else
        if sum(psi(:)) < 1
            psi = psi/sum(psi);
        end
    end
    
    % Relabel nodes from 1 to N
    E = relabel(E);
    
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
    
    msg  = cell(N,1);
    b    = cell(N,1);
    
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
    for i = 1:N
        tmp    = cell(size(nghb{i}));
        tmp(:) = {ones(card(i),1)};
        msg{i} = tmp;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  1st Phase: Leaves -> Root  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for idx = size(queue,1):-1:1
        orig = queue(idx,1);
        targ = queue(idx,2);
        
        blf = phi{orig}.*prod(cell2mat(msg{orig}), 2);
        
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
        
        msg{targ}{nghb2idx{targ}(orig)} = psi_cur'*(blf./msg{orig}{nghb2idx{orig}(targ)});  % here, parentheses DO matter! trust me!
        msg{targ}{nghb2idx{targ}(orig)}(isnan(msg{targ}{nghb2idx{targ}(orig)})) = 0;
        
        % Scale up if it is too small
        if sum(msg{targ}{nghb2idx{targ}(orig)}) < 1
            msg{targ}{nghb2idx{targ}(orig)} = msg{targ}{nghb2idx{targ}(orig)}./sum(msg{targ}{nghb2idx{targ}(orig)});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  2nd Phase: Root -> Leaves  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for idx = 1:size(queue,1)
        orig = queue(idx,2);
        targ = queue(idx,1);
                
        blf = phi{orig}.*prod(cell2mat(msg{orig}), 2);
        
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
        
        msg{targ}{nghb2idx{targ}(orig)} = psi_cur'*(blf./msg{orig}{nghb2idx{orig}(targ)});  % here, parentheses DO matter! trust me!
        msg{targ}{nghb2idx{targ}(orig)}(isnan(msg{targ}{nghb2idx{targ}(orig)})) = 0;
        
        % Scale up if it is too small
        if sum(msg{targ}{nghb2idx{targ}(orig)}) < 1
            msg{targ}{nghb2idx{targ}(orig)} = msg{targ}{nghb2idx{targ}(orig)}./sum(msg{targ}{nghb2idx{targ}(orig)});
        end        
    end
    
    %%%%%%%%%%%%%%%%%
    %%%  Beliefs  %%%
    %%%%%%%%%%%%%%%%%
    for i = 1:N
        b{i} = phi{i}.*prod(cell2mat(msg{i}), 2);
        b{i} = b{i}./sum(b{i});
    end
end