function comp = conn_comp(nghb)
% function comp = conn_comp(nghb)
% 
% This function returns the connected components of a graph defined by the
% list of neighbors nghb.
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
% 
% 
% EXAMPLES
%
% E.g. assume the simple graphical model:
%
%       1
%     /   \
%    2     3        5 -- 6
%           \
%            4
% 
% The nghb cell array is as follows:
% nghb{1} = [2 3];  nghb{2} = [1];  nghb{3} = [1 4];  nghb{4} = [3];
% nghb{5} = [6];    nghb{6} = [5];
%
% comp = conn_comp(nghb);
% 
% This will produce output:
% comp{1} = [1  2  3  4];  comp{2} = [5  6];
% 
% Author: geopapa
% $ Date: 2013/10/16 09:57:12 $

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
end