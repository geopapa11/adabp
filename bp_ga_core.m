function [mu_marg, S_marg, h_marg, J_marg, msg_h, msg_J, nghb] = bp_ga_core(h, J, varargin)
% function [mu_marg, S_marg, h_marg, J_marg, msg_h, msg_J, nghb] = bp_ga_core(h, J, varargin)
% 
% This function performs belief propagation in a gaussian tree mrf.
%
% INPUT
%        N:  1x1   scalar        - Number of nodes (variables)
%        d:  Nx1   double array  - Dimension of variables (d(i) is dimension for node i)
%        h:  Nx1 cell array      - Information vector
%            (if d > 1)            Each element is the information vector of
%            Nx1 array             the corresponding variable.
%            (if d = 1)
%        J:  NxN cell matrix     - Information (precision) matrix
%            (if d > 1)
%            NxN matrix
%            (if d = 1)
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
%  mu_marg:  Nx1 cell array      - Mean of the marginal
%                                  Each element is: x_i ~ N(mu_marg{i}, S_marg{i})
%   S_marg:  Nx1 cell array      - (Co)variance of the marginal
%                                  Each element is: x_i ~ N(mu_marg{i}, S_marg{i})
%   h_marg:  Nx1 cell array      - Information vector of the marginal
%                                  Each element is: x_i ~ N^{-1}(h_marg{i}, J_marg{i})
%   J_marg:  Nx1 cell array      - Information matrix of the marginal
%                                  Each element is: x_i ~ N^{-1}(h_marg{i}, J_marg{i})
%    msg_h:  Nx1 cell matrix     - Information messages
%                                  Each row i is a cell array representing
%                                  the information messages that are incoming to node i
%    msg_J:  Nx1 cell matrix     - Precision messages
%                                  Each row i is a cell array representing 
%                                  the precision messages that are incoming to node i
%     nghb:  Nx1 or 1xN          - Neighbors
%            cell array            Each element lists the neighbors of the
%                                  corresponding variable
% 
% Remember:  h = inv(Sigma)*mu,  J = inv(Sigma)
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
% with information vector (h) and information (precision) matrix (J):
% h = [1 1 1 1]';
%
% J = [1    0.8  0.2  0  ;
%      0.8  1    0    0  ;
%      0.2  0    1    0.3;
%      0    0    0.3  1  ];
%
% [mu_marg, S_marg, h_marg, J_marg] = bp_ga_core(h, J);
% 
% This will produce output:
% mu_marg = {[0.1460], [0.8832], [0.7371], [0.7789]}';
% S_marg  = {[3.1641], [3.0250], [1.2517], [1.1127]}';
% 
% Lastly, we make an example where different nodes (N=4) have different
% dimensions, d = [2 3 1 2].
% 
% h = {[1;1];  [1;1;1];  [1];  [1;1]};
% J = {[1 .3; .3 1],            [.1 .2 .1; .2 .1 -.1],          [.2;.1],   [0 0; 0 0];
% 	   [.1 .2; .2 .1; .1 -.1],  [1 -.2 .1; -.2 1 .3; .1 .3 1],  [0;0;0],   [0 0; 0 0; 0 0];
% 	   [.2 .1],                 [0 0 0],                        [1],       [.3 -.1];
% 	   [0 0; 0 0],              [0 0 0; 0 0 0],                 [.3; -.1], [1 -.4; -.4 1]};
%
% [mu_marg, S_marg, h_marg, J_marg] = bp_ga_core(h, J);
% 
% This will produce output:
% mu_marg{1} = [0.3633   0.6273]';  mu_marg{2} = [0.9403  0.8494  0.6776]';
% mu_marg{3} = [0.5481]';           mu_marg{4} = [1.4859  1.6528]';
% 
% S_marg{1} = [1.1950   -0.3051;  -0.3051   1.1926];
% S_marg{2} = [1.1557    0.3473   -0.2374;  0.3473    1.2504   -0.4077; -0.2374   -0.4077    1.1771];
% S_marg{3} = [1.1569];
% S_marg{4} = [1.3013    0.4847;   0.4847   1.1911];
% 
% Author: geopapa
% $ Date: 2013/10/16 16:31:03 $

    % Check if there is only one node in the graph
    if length(h)==1
        if iscell(h),  h = cell2mat(h);  end
        if iscell(J),  J = cell2mat(J);  end
        h_marg{1}  = h;
        J_marg{1}  = J;
        S_marg{1}  = inv(J_marg{1});
        mu_marg{1} = S_marg{1}*h_marg{1};
        msg_h      = {{zeros(size(h))}};
        msg_J      = {{zeros(size(J))}};
        return;
    end
    
    N    = length(h);
    nghb = cell(N,1);
    
    if (iscell(h) && ~iscell(J)) || (~iscell(h) && iscell(J))
        error('h and J have to be both either numeric or cell arrays.');
    end
    
    if ~iscell(h),  h = num2cell(h);  J = num2cell(J);  end
    
    d = num2cell(cellfun(@length,h));
    
    h_marg  = cell(N,1);
    J_marg  = cell(N,1);
    mu_marg = cell(N,1);
    S_marg  = cell(N,1);
    
    % Find the neighbors of each node
    if any(strcmp(varargin,'nghb'))
        nghb = varargin{find(strcmp(varargin,'nghb'))+1};
        nghb = relabel(nghb);
    else
        J_abs = cellfun(@abs,                  J, 'UniformOutput', false);
        idx   = cellfun(@(v)lt(v(:),1e-9), J_abs, 'UniformOutput', false);
        idx   = ~cellfun(@all,idx);
        idx   = logical(idx - eye(size(idx)));
        for i = 1:N
            nghb{i} = find(idx(i,:));
        end
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
    
    % Generate the schedule of messages to be sent
    queue = msg_schedule(root, ch);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Initialize the messages  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg_h = cell(N,1);
    msg_J = cell(N,1);
    for i = 1:N
        tmp      = cell(size(nghb{i}));
        tmp(:)   = {zeros(d{i},1)};
        msg_h{i} = tmp;
        tmp(:)   = {zeros(d{i})};
        msg_J{i} = tmp;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  1st Phase: Leaves -> Root  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for idx = size(queue,1):-1:1
        orig = queue(idx,1);
        targ = queue(idx,2);
        
        h_marg_tmp = h{orig}      + sum(cat(3,msg_h{orig}{:}),3);
        J_marg_tmp = J{orig,orig} + sum(cat(3,msg_J{orig}{:}),3);
        
        S = -J{targ,orig}/(J_marg_tmp - msg_J{orig}{nghb2idx{orig}(targ)});
        msg_h{targ}{nghb2idx{targ}(orig)} = S*(h_marg_tmp - msg_h{orig}{nghb2idx{orig}(targ)});
        msg_J{targ}{nghb2idx{targ}(orig)} = S*J{orig,targ};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  2nd Phase: Root -> Leaves  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for idx = 1:size(queue,1)
        orig = queue(idx,2);
        targ = queue(idx,1);
        
        h_marg_tmp = h{orig}      + sum(cat(3,msg_h{orig}{:}),3);
        J_marg_tmp = J{orig,orig} + sum(cat(3,msg_J{orig}{:}),3);
        
        S = -J{targ,orig}/(J_marg_tmp - msg_J{orig}{nghb2idx{orig}(targ)});
        msg_h{targ}{nghb2idx{targ}(orig)} = S*(h_marg_tmp - msg_h{orig}{nghb2idx{orig}(targ)});
        msg_J{targ}{nghb2idx{targ}(orig)} = S*J{orig,targ};
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Beliefs in Information Form  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N
        h_marg{i} = h{i}   + sum(cat(3,msg_h{i}{:}),3);
        J_marg{i} = J{i,i} + sum(cat(3,msg_J{i}{:}),3);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Beliefs in Moment Form  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N
        S_marg{i}  = inv(J_marg{i});
        mu_marg{i} = S_marg{i}*h_marg{i};
    end
end