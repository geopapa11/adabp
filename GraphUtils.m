classdef GraphUtils < handle
    %INCRBP Summary of this class goes here
    %   Detailed explanation goes here
    % Author: geopapa
    % $ Date: 2014/02/15 11:12:22 $

    
   methods (Static = true)
       %% This function relabels the sets to new labels specified by argument lbl.
       %  If no argument is specified, the new labeling is sequential, in other
       %  words, new labeling starts from label 1 and ends in label N, where N 
       %  is the number of distinct nodes (elements) of the sets.
       %
       % INPUT
       %        N:  1x1   scalar      - Number of distinct nodes (variables)
       %        C:  1x1   scalar      - Number of cliques
       %     sets:  Cx1 cell array    - Collection of sets
       %            or                  Each cell element is a set
       %            CxM double array    Each element is a set of size M
       %                                
       %      lbl:  Nx1    array      - Labels (the new labels)
       %      (optional)                
       % 
       % 
       % OUTPUT
       %     sets:  Cx1 cell array    - Collection of sets
       %            or                  Each cell element is a set
       %            CxM double array    Each element is a set of size M
       %                                (Nodes are relabeled as specified by lbl
       %                                or from 1 to N, where N is the number of
       %                                distinct nodes)       
       function sets_out = relabel(sets, lbl)
           sets_out = sets;
           if iscell(sets_out)
               keys = unique(cell2mat(sets_out'));
           else
               keys = unique(sets_out(:));
           end
           
           if nargin == 2
               vals = unique(lbl,'stable');
           else
               vals = 1:length(keys);
           end
           
           if length(keys) ~= length(vals)
               error('The new labels must be as many as the old labels.');
           end
           
           map = sparse(keys,1,vals,max(keys),1);
           
           if iscell(sets_out)
               for k = 1:length(sets_out)
                   sets_out{k} = full(map(sets_out{k}));
                   if iscolumn(sets_out{k}),  sets_out{k} = sets_out{k}';  end
               end
           else
               sets_out = full(map(sets_out));
               if iscolumn(sets_out),  sets_out = sets_out';  end
           end
       end  % end relabel
       
       function nghb = find_nghb(E)
           if isempty(E)
               nghb = [];               
           else
               E    = relabel(E);
               N    = length(unique(E(:)));
               nghb = cell(N,1);
               for i = 1:N
                   [row,col]   = find(E==i);
                   col         = col+1;
                   col(col==3) = 1;
                   idx         = sub2ind(size(E),row,col);
                   nghb{i}     = E(idx)';
               end
               
               % Sort neighbors w.r.t. their indices
               nghb = cellfun(@sort, nghb, 'UniformOutput', false);
           end
       end  % end find_nghb
       
       %% This function returns the children of each node based on the chosen root
       %  and the list of neighbors nghb.
       %
       % INPUT
       %      N:  1x1   scalar      -  Number of nodes (variables)
       %   nghb:  Nx1 or 1xN        -  Neighbors
       %          cell array           Each element lists the neighbors of the
       %                               corresponding variable.
       %   root:  1x1   scalar      -  Root
       % 
       % 
       % OUTPUT
       %     ch:  Nx1 cell array    -  Children
       %                               Each element, which is a 1x(num_of_children)
       %                               array, lists the children of each node based
       %                               on the chosen root.
       function ch = get_ch(nghb,root)
           N  = length(nghb);
           
           % Choose a root randomly if not specified
           if nargin == 1,  root = randi(N);  end
                      
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
       end  % end get_ch
       
       %% This function returns the Eulerian trail of a graph, in other words, 
       %  a trail which visits every edge exactly once. 
       %
       % INPUT
       %   root:  1x1   scalar        -  Root
       %     ch:  Nx1 cell array      -  Children
       %                                 Each element, which is a 1x(num_of_children)
       %                                 array, lists the children of each node based
       %                                 on the chosen root.
       % 
       % OUTPUT
       %     E:  1x2N-1 double array  -  Euler tour
       %    De:  1x2N-1 double array  -  Depth of each node in the Euler tour
       %    pa:  1xN    double array  -  Parent of a node
       %     D:  1xN    double array  -  Depth of a node
       %     H:  1xN    double array  -  Index of the first occurence of node i in E
       function [E, De, H, D, pa] = eulertour(root, ch)
           N   = length(ch);
           E   = zeros(1,2*N-1);  % Euler Tour
           De  = zeros(1,2*N-1);  % Depth of each node in the Euler Tour
           pa  = zeros(1,N);      % Parents of each node
           D   = zeros(1,N);      % Depth of each node
           H   = zeros(1,N);      % Index of the first occurence of node i in E
           dfs = zeros(1,N);      % DFS sequence path
           
           % Traverse the tree in a depth-first manner, and store the parent and
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
       end  % end eulertour
       
       %% This function returns the connected components of a graph defined by the
       %  list of neighbors nghb.
       %
       % INPUT
       %      N:  1x1   scalar      - Number of nodes (variables)
       %      M:  1x1   scalar      - Number of connected components
       %   nghb:  Nx1 cell array    - Neighbors
       %                              Each element lists the neighbors of the
       %                              corresponding variable.
       % 
       % 
       % OUTPUT
       %     comp:  Mx1 cell array  - Connected components
       %                              Each element lists the nodes that constitute
       %                              the corresponding connected component.
       %                              (Remember, a connected component is a
       %                              subgraph where every pair of its nodes are
       %                              reachable through some path)       
       function comp = conn_comp(nghb)
           N         = length(nghb);
           cnt       = 1;
           stack     = 1;
           expl_set  = 1:N;
           comp{cnt} = 1;
           while true        
               node       = stack(end);
               stack(end) = [];
               nghbs      = nghb{node};
               
               if ~isempty(nghbs)
                   for k = nghbs
                       idx          = nghb{k}==node;
                       nghb{k}(idx) = [];
                   end
                   
                   comp{cnt} = [comp{cnt}  nghbs];
                   stack     = [stack      nghbs];
                   
                   for k = nghbs
                       expl_set(expl_set==k) = [];
                   end
               end
               
               expl_set(expl_set==node) = [];
               
               if isempty(expl_set)
                   break;
               end
               
               if isempty(stack)
                   cnt       = cnt + 1;
                   comp{cnt} = expl_set(1);
                   stack     = expl_set(1);
               end
           end
           
           comp = cellfun(@unique, comp, 'UniformOutput', false);
           comp = comp';
       end  % end conn_comp
   end % methods
end % classdef