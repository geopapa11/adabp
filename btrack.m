function seq = btrack(root, root_idx, nghb, msg_idx, varargin)
% function seq = btrack(root, root_idx, nghb, msg_idx, varargin)
% 
% This function recovers the maximum sequence of values by backtracking.
%
% INPUT
%         N:  1x1   scalar     -  Number of nodes (variables)
%      root:  1x1   scalar     -  Root
% root_idx:  1x1   scalar      -  Index corresponding to the value of the
%                                 chosen root
%      nghb:  Nx1 cell array   -  Neighbors
%                                 Each element lists the neighbors of the
%                                 corresponding variable
%                                 If the N nodes are not labeled from 1 to
%                                 N (and have arbitrary labels), then the
%                                 first element in cell array nghb is
%                                 expected to correspond to the node with
%                                 the smallest label and the following
%                                 elements (in cell array nghb) to nodes
%                                 with increasing labels.
%  msg_idx:  Nx1 cell array    -  Message maximum indices
%                                 Each row lists the indices of the source
%                                 nodes that maximize the value of the
%                                 message at the target (i) node
%  varargin: cell array        -  Contains pairs of arguments
%                                 * If varargin{i} = 'nghb2idx', then varargin{i+1}
%                                   is a containers.Map representing the map
%                                   of neighbors of a node to their linear
%                                   indices
%                                 * If varargin{i} = 'seq_prv', then varargin{i+1}
%                                   is the sequence obtained at the
%                                   previous iteration
%                                 * If varargin{i} = 'replace_nghb', then varargin{i+1}
%                                   is a 1x2 cell array. 
%                                   First element is the index of the node
%                                   whose neighbors will be set to a new
%                                   value and the second element is the new
%                                   neighbors of this node.
%                                 * If varargin{i} = 'isDirty', then varargin{i+1}
%                                   is a NxN sparse logical matrix. 
%                                   Row indicates the source and column the
%                                   target node.
%                                   If isDirty(i,j)=1, this means that the 
%                                   delta message i -> j has changed from
%                                   the previous iteration.
% 
% 
% OUTPUT
%      seq:  1xN double array  -  Sequence of variable (node) values
%                                 corresponding to the optimal value
% 
% 
% EXAMPLES
%
% E.g., assume the simple graphical model:
%
%       1
%     /   \
%    2     3
%           \
%            4
% 
% with root, map_root, and neighbors and message maximum indices, as
% follows:
% root = 1;
% root_idx = 1;
% 
% nghb{1} = [2  3];  nghb{2} = [1];  nghb{3} = [1  4];  nghb{4} = [3];
% 
% msg_idx    = cell(4,1);
% msg_idx{1} = {[2;2], [1;1]};
% msg_idx{2} = {[2;1;1]};
% msg_idx{3} = {[1;1], [3;2]};
% msg_idx{4} = {[1;2;1;2]};
% 
% seq = btrack(root, root_idx, nghb, msg_idx);
% 
% This will produce output:
% seq = [1  2  1  3];
% 
% Now, let's assume that nodes in the graph have arbitrary labeling. (This
% can arise in situations where we consider a connected component from a
% forest):
% 
%       8
%     /   \
%    4     9
%           \
%            7
% 
% with root, map_root, and neighbors and message maximum indices, as
% follows:
% root = 8;
% root_idx = 1;
% 
% nghb{1} = [8];  nghb{2} = [9];  nghb{3} = [4  9];  nghb{4} = [8  7];
% 
% msg_idx    = cell(4,1);
% msg_idx{1} = {[2;1;1]};
% msg_idx{2} = {[1;2;1;2]};
% msg_idx{3} = {[2;2], [1;1]};
% msg_idx{4} = {[1;1], [3;2]};
% 
% seq = btrack(root, root_idx, nghb, msg_idx);
% 
% This will produce output:
% seq = [2  3  1  1];
% 
% Here, the first elements of nghb and msg_idx correspond to the node with
% the smallest labeling, that is, "4". The second ones to the node with the
% second smallest labeling, that is, "7". The third ones to the node with
% label "8" and the last ones to the node with the largest label, "9".
% 
% Author: geopapa
% $ Date: 2014/05/17 12:56:48 $

    N = length(nghb);

    % Create a map from node labels to linear indices 1:N, where N is the # of nodes
    if any(strcmp(varargin,'lbl2idx'))
        lbl2idx = varargin{find(strcmp(varargin,'lbl2idx'))+1};
    else
        lbl     = unique([nghb{:}]);
        lbl2idx = sparse(lbl,1,1:N,max(lbl),1);
    end
    
    if nnz(lbl2idx) ~= N
        save('debugging.mat', 'root', 'root_idx', 'nghb');
        error('The number of distinct labels should equal the number of distinct nodes.');
    end
    
    % Create a map of neighbors to their linear indices in cell matrix nghb
    if any(strcmp(varargin,'nghb2idx'))
        nghb2idx = varargin{find(strcmp(varargin,'nghb2idx'))+1};
    else
        nghb2idx = cell(N,1);
        for i = 1:N
            keys        = nghb{i};
            vals        = 1:length(nghb{i});
            nghb2idx{i} = sparse(keys,1,vals,max(keys),1);
        end
    end
    
    % If the sequence at the previous iteration is known and the maximizing
    % value of a node stayed the same, then the values of all the nodes of
    % the subtree rooted at this node will remain the same (provided the
    % delta messages did not change). As a consequence, there is no need to
    % re-compute these values.
    if any(strcmp(varargin,'seq_prv'))
        seq_prv      = varargin{find(strcmp(varargin,'seq_prv'))+1};
        seq          = seq_prv;
        prvSeqExists = true;
    else
        seq          = zeros(1,N);
        prvSeqExists = false;
    end
    
    % This argument provides the flexibility to change the neighbors of a
    % single node
    if any(strcmp(varargin,'replace_nghb'))
        tmp       = varargin{find(strcmp(varargin,'replace_nghb'))+1};
        idx       = tmp{1};
        nghb_upd  = tmp{2};
        nghb{idx} = nghb_upd;
    end
    
    if prvSeqExists
        if any(strcmp(varargin,'isDirty'))
            isDirty = varargin{find(strcmp(varargin,'isDirty'))+1};
        else
            isDirty = logical(sparse(N,N));
        end
    end
    
    if prvSeqExists
        seq = assign_prvSeq(seq, seq_prv, root, root_idx, lbl2idx, nghb, nghb2idx, msg_idx, isDirty);
    else
        seq = assign(seq, root, root_idx, lbl2idx, nghb, nghb2idx, msg_idx);
    end
end

function seq = assign(seq, root, root_idx, lbl2idx, nghb, nghb2idx, msg_idx)
    i = root;
    seq(lbl2idx(i)) = root_idx;
    stack = i;
    while ~isempty(stack)
        i          = stack(end);
        stack(end) = [];
        for j = nghb{lbl2idx(i)}
            seq(lbl2idx(j)) = msg_idx{lbl2idx(i)}{nghb2idx{lbl2idx(i)}(j)}(seq(lbl2idx(i)));  % Maximizing indices of x_j for every value of x_i
            nghb{lbl2idx(j)}(nghb{lbl2idx(j)}==i) = [];
        end
        stack = [stack, nghb{lbl2idx(i)}];
    end
end

function seq = assign_prvSeq(seq, seq_prv, root, root_idx, lbl2idx, nghb, nghb2idx, msg_idx, isDirty)
    i = root;
    seq(lbl2idx(i)) = root_idx;
    stack = i;
    while ~isempty(stack)
        i          = stack(end);
        stack(end) = [];
        for j = nghb{lbl2idx(i)}
            if ( seq(lbl2idx(i)) ~= seq_prv(lbl2idx(i)) ) || isDirty(j,i)
                seq(lbl2idx(j)) = msg_idx{lbl2idx(i)}{nghb2idx{lbl2idx(i)}(j)}(seq(lbl2idx(i)));  % Maximizing indices of x_j for every value of x_i
                nghb{lbl2idx(j)}(nghb{lbl2idx(j)}==i) = [];
                stack = [stack, j];
            end
        end
    end
end