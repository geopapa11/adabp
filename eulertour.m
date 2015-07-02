function [E, De, H, D, pa] = eulertour(root, ch)
% root   = 1;
% 
% ch     = cell(1,10);
% ch{1}  = [2 3];
% ch{2}  = [4 5];
% ch{3}  = [6 7];
% ch{4}  = [];
% ch{5}  = [8];
% ch{6}  = [];
% ch{7}  = [9 10];
% ch{8}  = [];
% ch{9}  = [];
% ch{10} = [];

% Euler tour and depth
    N = length(ch);
    
    E  = zeros(1,2*N-1);  % Euler Tour
    De = zeros(1,2*N-1);  % Depth of each node in the Euler Tour    
    
    pa  = zeros(1,N);     % Parents of each node
    D   = zeros(1,N);     % Depth of each node
    H   = zeros(1,N);     % Index of the first occurence of node i in E
    dfs = zeros(1,N);     % DFS sequence path
    
    % Traverse the tree in a depth-first-manner, and store the parent and
    % depth of each node
    
    stack = java.util.Stack;
    stack.push(root);
    cnt = 1;
    while ~stack.isEmpty
        node = stack.pop;        
        % Push the children onto the stack
        for i = length(ch{node}):-1:1
            stack.push(ch{node}(i));
        end
        
        % Determine the depth of nodes
        for i = 1:length(ch{node})
            pa(ch{node}) = node;
            D(ch{node})  = D(node)+1;
        end
        
        % Store the current node
        dfs(cnt) = node;
        cnt      = cnt + 1;        
    end
    
    E(1) = dfs(1);
    cnt  = 2;
    for j = 2:N   % index of the dfs path
        prev = dfs(j-1);        
        while pa(dfs(j)) ~= prev
            prev   = pa(prev);
            E(cnt) = prev;
            cnt    = cnt + 1;
        end        
        E(cnt) = dfs(j);
        cnt    = cnt + 1;
    end
    
    if N > 1
        prev = pa(dfs(N)); 
        while prev ~= root
            E(cnt) = prev;
            cnt    = cnt + 1;
            prev   = pa(prev);
        end
        E(cnt) = root;
    end
    
    for k = 1:2*N-1  % index of the euler tour
        De(k) = D(E(k));
    end
    
    for k = 1:2*N-1
        if H(E(k)) == 0
            H(E(k)) = k;
        end
    end
end