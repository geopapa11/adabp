function queue = msg_schedule(root, ch)
% function queue = msg_schedule(root, ch)
% 
% This function creates the message schedule starting from the root and
% propagating down to the leaves.
%
% INPUT
%      N:  1x1   scalar      - Number of nodes (variables)
%   root:  1x1   scalar      - Root
%     ch:  Nx1 cell array    - Children
%                              Each element lists the children of each node 
%                              based on the chosen root.
% 
% 
% OUTPUT
%  queue:  (N-1)x2 matrix    - Message queue
%                              Each row has two elements. The first element
%                              denotes the origin node and the second the
%                              target node.
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
% with root:
% root = 1;
%
% The ch cell array will be:
% ch{1} = [2 3];  ch{2} = zeros(1,0);  ch{3} = [4];  ch{4} = zeros(1,0);
% 
% queue = msg_schedule(root, ch);
% 
% This will produce output:
% queue = [3  1;
%          2  1;
%          4  3]
% 
% That is, we start from the messages that originate from the root (1) and
% propagate downwards.
% 
% Author: geopapa
    
    for i = 1:length(ch)
        ch{i} = ch{i}(:)';
    end

    N     = length(ch);    
    i     = root;
    count = 1;
    stack = [];
    queue = zeros(N-1,2);
    while count <= N-1
        for k = length(ch{i}):-1:1
            queue(count, :) = [ch{i}(k), i];
            count           = count + 1;
        end        
        stack      = [stack, ch{i}];
        i          = stack(end);
        stack(end) = [];
    end    
end
