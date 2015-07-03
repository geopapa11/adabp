function sets_out = relabel(sets, lbl)
% function sets = relabel(sets)
% 
% This function relabels the sets to new labels specified by argument lbl.
% If no argument is specified, the new labeling is sequential, in other
% words, new labeling starts from label 1 and ends in label N, where N 
% is the number of distinct nodes (elements) of the sets.
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
% 
% 
% EXAMPLES
%
% E.g. assume the simple factor graph:
%
%    2    4   8
%    | \  |  /|
%    |  \ | / |
%    c    a   d
% 
% with sets:
% sets{1} = [2 8 4]';
% sets{2} = [2]';
% sets{3} = [8]';
% 
% sets = relabel(sets);
% 
% This will produce output:
% sets{1} = [1 3 2]';
% sets{2} = [1]';
% sets{3} = [3]';
% 
% 
% This is an example of cliques of fixed size M (here, M=2).
% 
%       3
%     /   \
%    4     7
%           \
%            9
% 
% The sets are:
% sets = [3  4;
%         3  7;
%         7  9];
% 
% sets = relabel(sets);
% 
% This will produce output:
% sets = [1  2;
%         1  3;
%         3  4];
% 
% Lastly, if you have:
%       1
%     /   \
%    2     3
%           \
%            4
% 
% and you want to relabel it as specified by lbl = [3 4 7 9], this will
% produce the following output:
% sets = [1  2;
%         1  3;
%         3  4];
% 
% sets = relabel(sets,lbl);
% 
% sets = [3  4;
%         3  7;
%         7  9];
% 
% Author: geopapa
% $ Date: 2013/02/14 18:21:33 $

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
end