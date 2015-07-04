function [mu_marg, S_marg, h_marg, J_marg, msg_h, msg_J, nghb] = bp_ga(h, J, obs, varargin)
% function [mu_marg, S_marg, h_marg, J_marg, msg_h, msg_J, nghb] = bp_ga(h, J, obs, varargin)
% 
% This function performs belief propagation in a gaussian mrf with 
% connected comonents and possibly observed variables.
%
% INPUT
%        N:  1x1   scalar        -  Number of nodes (variables)
%        d:  Nx1   double array  -  Dimension of variables (d(i) is dimension for node i)
%        h:  Nx1 cell array      -  Information vector
%            (if d > 1)             Each element is the information vector of
%            Nx1 array              the corresponding variable.
%            (if d = 1)
%        J:  NxN cell matrix     -  Information (precision) matrix
%            (if d > 1)          
%            NxN matrix          
%            (if d = 1)
%      obs:  Nx1 cell array      -  Observed values
%            (if d > 1)             Each element equals to the value that the
%            Nx1 array              corresponding variable is set to.
%            (if d = 1)             (The element(s) of an unobserved variable 
%                                   is (are) set to NaN)
%  varargin: cell array          -  Contains pairs of arguments
%                                    * If varargin{i} = 'root', then varargin{i+1}
%                                      is an 1x1 scalar representing the root.
%                                      (Default value: root chosen randomly)
% 
% 
% OUTPUT
%  mu_marg:  Nx1 cell array    -  Mean of the marginal
%                                 Each element is: x_i ~ N(mu_marg{i}, S_marg{i})
%   S_marg:  Nx1 cell array    -  (Co)variance of the marginal
%                                 Each element is: x_i ~ N(mu_marg{i}, S_marg{i})
%   h_marg:  Nx1 cell array    -  Information vector of the marginal
%                                 Each element is: x_i ~ N^{-1}(h_marg{i}, J_marg{i})
%   J_marg:  Nx1 cell array    -  Information matrix of the marginal
%                                 Each element is: x_i ~ N^{-1}(h_marg{i}, J_marg{i})
%    msg_h:  Nx1 cell matrix   -  Information messages
%                                 Each row i is a cell array representing
%                                 the information messages that are incoming to node i
%    msg_J:  Nx1 cell matrix   -  Precision messages
%                                 Each row i is a cell array representing 
%                                 the precision messages that are incoming to node i
%     nghb:  Nx1 or 1xN        -  Neighbors
%            cell array           Each element lists the neighbors of the
%                                 corresponding variable
% 
% Remember:  h = inv(Sigma)*mu,  J = inv(Sigma)
% 
% 
% EXAMPLES
%
% E.g. assume graphical model comprised of two connected tree components:
% 
%       1
%     /   \
%    2     3        5 -- 6
%           \
%            4
% 
% with information vector and information (precision) matrix as follows:
% h = [1 1 1 1 1 1]';
%
% J = [1    0.8  0.2  0    0    0  ;
%      0.8  1    0    0    0    0  ;
%      0.2  0    1    0.3  0    0  ;
%      0    0    0.3  1    0    0  ;
%      0    0    0    0    1    0.7;
%      0    0    0    0    0.7  2  ];
%
% [mu_marg, S_marg, h_marg, J_marg] = bp_ga(h, J);
% 
% This will produce output:
% 
% mu_marg = {[0.1460], [0.8832], [0.7371], [0.7789], [0.8609], [0.1987]}';
% S_marg  = {[3.1641], [3.0250], [1.2517], [1.1127], [1.3245], [0.6623]}';
% 
% 
% Lastly, we show an example with observed variables.
% 
%       1
%     /   \
%    2     3        5 -- 6
%           \
%            4
% 
% with information vector, information (precision) matrix and observations as follows:
% h = [1 1 1 1 1 1]';
%
% J = [1    0.8  0.2  0    0    0  ;
%      0.8  1    0    0    0    0  ;
%      0.2  0    1    0.3  0    0  ;
%      0    0    0.3  1    0    0  ;
%      0    0    0    0    1    0.7;
%      0    0    0    0    0.7  2  ];
% 
% obs = [NaN  NaN  4  NaN NaN -1]';
%
% [mu_marg, S_marg, h_marg, J_marg] = bp_ga(h, J, obs);
% 
% This will produce output:
% 
% mu_marg = {[-1.6667], [2.3333], [ 4], [-0.2000], [1.7000], [ -1]}';
% S_marg  = {[2.7778] , [2.7778], [ 0], [ 1]     , [ 1]    , [  0]}';
% 
% Author: geopapa
% $ Date: 2013/10/16 18:43:28 $
    
    if (iscell(h) && ~iscell(J)) || (~iscell(h) && iscell(J))
        error('h and J have to be both either numeric or cell arrays.');
    end
    
    if ~iscell(h),  h = num2cell(h);  J = num2cell(J);  end
    
    d = num2cell(cellfun(@length,h));
    
    N       = length(h);
    h_marg  = cell(N,1);
    J_marg  = cell(N,1);
    mu_marg = cell(N,1);
    S_marg  = cell(N,1);
    msg_h   = cell(N,1);
    msg_J   = cell(N,1);
    
    % Check if there are observed variables
    if nargin >= 3 && ~isempty(obs)
        if ~iscell(obs)
            obs = num2cell(obs);
        end
        
        any_nan_idx = cellfun(@(v) any(isnan(v(:))), obs);
        all_nan_idx = cellfun(@(v) all(isnan(v(:))), obs);
        if any(any_nan_idx ~= all_nan_idx)
            error('If a variable is observed, we make the assumption that all its dimensions should be observed.');
        end
        
        unobs_idx      = all_nan_idx;
        obs(unobs_idx) = [];
        
        h_cnd = cell2mat(h(unobs_idx)) - cell2mat(J(unobs_idx,~unobs_idx))*cell2mat(obs);
        h_cnd = mat2cell(h_cnd, cell2mat(d(unobs_idx)));
        J_cnd = J(unobs_idx,unobs_idx);
        h     = h_cnd;
        J     = J_cnd;
    else
        unobs_idx = true(1,N);
    end
    
    % Determine the root
    if ~any(strcmp(varargin,'root'))
        root = randi(sum(unobs_idx));  % Choose a root randomly from the unobserved nodes
    else
        root          = varargin{find(strcmp(varargin,'root'))+1};
        unobs_lin_idx = find(unobs_idx);
        root_unobs    = find(unobs_lin_idx==root, 1);
        if isempty(root_unobs),  error('You chose a node as a root that is observed.');  end
        root = root_unobs;
    end
    
    N_uno       = length(h);      % number of unobserved variables
    nghb        = cell(N_uno,1);
    h_marg_uno  = cell(N_uno,1);
    J_marg_uno  = cell(N_uno,1);
    mu_marg_uno = cell(N_uno,1);
    S_marg_uno  = cell(N_uno,1);
    msg_h_uno   = cell(N_uno,1);
    msg_J_uno   = cell(N_uno,1);
    
    % Find the neighbors of each node
    J_abs = cellfun(@abs,                  J, 'UniformOutput', false);
    idx   = cellfun(@(v)lt(v(:),1e-9), J_abs, 'UniformOutput', false);
    idx   = ~cellfun(@all,idx);
    idx   = logical(idx - eye(size(idx)));
    for i = 1:N_uno
        nghb{i} = find(idx(i,:));
    end
    
    % Get the connected components of the graph
    comp = conn_comp(nghb);
    
    % Run gaussian bp for each connected component
    if length(comp)==1
        [mu_marg_uno, S_marg_uno, h_marg_uno, J_marg_uno, msg_h_uno, msg_J_uno] = bp_ga_core(h, J, 'root', root, 'nghb', nghb);
    else
        for i = 1:length(comp)
            h_cc    = h(comp{i});
            J_cc    = J(comp{i},comp{i});
            nghb_cc = nghb(comp{i});
            
            [mu_marg_cc, S_marg_cc, h_marg_cc, J_marg_cc, msg_h_cc, msg_J_cc] = bp_ga_core(h_cc, J_cc, 'nghb', nghb_cc);
            
            mu_marg_uno(comp{i}) = mu_marg_cc;
            S_marg_uno(comp{i})  = S_marg_cc;
            h_marg_uno(comp{i})  = h_marg_cc;
            J_marg_uno(comp{i})  = J_marg_cc;
            msg_h_uno(comp{i})   = msg_h_cc;
            msg_J_uno(comp{i})   = msg_J_cc;
        end
    end
    
    if N_uno < N
        mu_marg(unobs_idx) = mu_marg_uno;
        S_marg(unobs_idx)  = S_marg_uno;
        h_marg(unobs_idx)  = h_marg_uno;
        J_marg(unobs_idx)  = J_marg_uno;
        msg_h(unobs_idx)   = msg_h_uno;
        msg_J(unobs_idx)   = msg_J_uno;
        
        mu_marg(~unobs_idx) = obs;
        S_marg(~unobs_idx)  = cellfun(@zeros,        d(~unobs_idx), 'UniformOutput', false);
        h_marg(~unobs_idx)  = cellfun(@(x) inf(x,1), d(~unobs_idx), 'UniformOutput', false);
        J_marg(~unobs_idx)  = cellfun(@inf,          d(~unobs_idx), 'UniformOutput', false);
        
        tmp            = cell(N,1);
        tmp(unobs_idx) = nghb;
        nghb           = tmp;
    else
        mu_marg = mu_marg_uno;
        S_marg  = S_marg_uno;
        h_marg  = h_marg_uno;
        J_marg  = J_marg_uno;
        msg_h   = msg_h_uno;
        msg_J   = msg_J_uno;
    end
end