function ch = get_ch(root, nghb)
% function ch = get_ch(root, nghb)
% 
% This function returns the children of each node based on the chosen root
% and the list of neighbors nghb.
%
% INPUT
%      N:  1x1   scalar      -  Number of nodes (variables)
%   root:  1x1   scalar      -  Root
%   nghb:  Nx1 or 1xN        -  Neighbors
%          cell array           Each element lists the neighbors of the
%                               corresponding variable.
% 
% 
% OUTPUT
%     ch:  Nx1 cell array    -  Children
%                               Each element, which is a 1x(num_of_children)
%                               array, lists the children of each node based
%                               on the chosen root.
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
% with root and neighbors as follows:
% root = 1;
%
% nghb{1} = [2 3];  nghb{2} = [1];  nghb{3} = [1 4];  nghb{4} = [3];
%
% ch = get_ch(root, nghb);
% 
% This will produce output:
% ch{1} = [2 3];
% ch{2} = [];
% ch{3} = [4];
% ch{4} = [];
% 
% Author: geopapa
% $ Date: 2013/10/15 14:12:04 $

    N     = length(nghb);
    ch    = cell(N,1);
    node  = root;
    cnt   = 1;
    stack = [];
    while cnt <= N
        ch{node} = nghb{node};
        for child = ch{node}
            idx = nghb{child}==node;
            nghb{child}(idx) = [];
        end
        stack = [stack, ch{node}];
        if isempty(stack)
            break;
        end
        node       = stack(end);
        stack(end) = [];
        cnt        = cnt + 1;
    end
end