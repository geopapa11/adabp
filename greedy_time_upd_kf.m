function [gr_set, f_gr, ig_gr, incr_gr, cnd_ig_gr, mean_up, cov_up, mean_gr, cov_gr, cov_wlk] = greedy_time_upd_kf(suffx, walk, set, expl, Y, c, budget, lambda, A, muV, Q, C, muW, R, mu1, cov1, flagCovWlk, varargin)
% function [gr_set, f_gr, ig_gr, incr_gr, cnd_ig_gr, mean_up, cov_up, mean_gr, cov_gr, cov_wlk] = greedy_time_upd_kf(suffx, walk, set, expl, Y, c, budget, lambda, A, muV, Q, C, muW, R, mu1, cov1, flagCovWlk, varargin)
%
% This function finds the greedy set, its reward, MI, the incremental reward
% and conditional MI at each step for the specified visit walk.
% 
% Lastly, it gives the greedy estimates of the parameter of interest:
% E(X)_{k|T}   = { E[X_1 | y_{1:T}], E[X_2 | y_{1:T}], ..., E[X_T | y_{1:T}] }
% and
% cov(X)_{k|T} = { cov(X_1 | y_{1:T}), cov(X_2 | y_{1:T}), ..., cov(X_T | y_{1:T}) }.
% 
% This function properly updates the initial mean and covariance of the 
% next element of the walk and then runs Kalman Filter only on the next 
% element of the walk (whose initial mean and covariance were properly 
% updated in the previous step).
% 
%
% INPUT
% 
%           d:  1x1  scalar        -  Dimension of state space
%           m:  1x1  scalar        -  Dimension of observed space
%           S:  1x1  scalar        -  Length of the visit walk
%           T:  1x1  scalar        -  Number of time points
%        walk:  1xS  array         -  Visit walk
%         set:  1xS  array         -  Measurements chosen for the corresponding visit walk
%                                     (if empty, measurement set will be chosen greedily)
%                                     (if an element is zero, this indicates absence of a measurement)
%        expl:  1xT  cell matrix   -  Exploration sets (for each time point)
%                                     If cell array or some of its elements are empty, 
%                                     then all measurements will be explored for the corresponding time points
%                                     (each element is an 1xN_j array denoting the measurement indices for each observation set)
%           Y:  1xT  cell matrix   -  Observed sequence (each element has up to m elements)
%           c:  mxT  matrix        -  Cost of each observation
%      budget:  1x1  scalar        -  Budget (positive number)
%      lambda:  1x1  scalar        -  Regularization parameter
%           A:  1xT cell matrix    -  Scaling of the state space at the previous time point
%                                     (each element is a dxd matrix)
%         muV:  1x(T-1) cell array -  Mean of noise process V at time k
%                                     (each element is a dx1 matrix)
%           Q:  1x(T-1) cell array -  Covariance of noise process V at time k
%                                     (each element is a dxd matrix)
%           C:  1xT cell array     -  Scaling factor of hidden states
%                                     (each element is a mxd matrix)
%         muW:  1xT cell array     -  Mean of noise process W at time k
%                                     (each element is an mx1 matrix)
%           R:  1xT cell array     -  Covariance of the noise term w_t (different for each t)
%                                     (each element is an mxm matrix)
%         mu1:  dx1  array         -  Mean of state_seq(1)
%        cov1:  dxd  matrix        -  Covariance of state_seq(1)
%  flagCovWlk:  1x1 boolean        -  Flag indicating whether the covariance
%                                     updates for each step of the walk would be generated
% 
% 
% OUTPUT
% 
%      gr_set:  1xS  array       -  Greedy set showing which measurements have been 
%                                   chosen greedily for each point in the walk
%        f_gr:  1x1  scalar      -  Reward of the greedy solution
%       ig_gr:  1x1  scalar      -  MI of the greedy solution
%     incr_gr:  1xT  array       -  Incremental reward of the greedy solution for each time point
%                                   incr_gr(k) = I(X; Y_k | y_{1:k-1}) - lambda*c_k
%  cnd_ig_gr:  1xT  array       -  Conditional MI of the greedy solution for each time point
%                                   cnd_ig_gr(k) = I(X; Y_k | y_{1:k-1})
%     mean_up:  dxT  matrix      -  E(X)_{k|k}    = {E[x_1 | y_{1}], E[x_2 | y_{1:2}], ..., E[x_T | y_{1:T}]}
%      cov_up:  1xT cell matrix  -  cov(X)_{k|k}  = cov(X_k | y_{1:k})
%                                   (each element is a dxd matrix)
%     mean_gr:  dxT  matrix      -  E(X)_{k|T}    = {E[x_1 | y_{1:T}], E[x_2 | y_{1:T}], ..., E[x_T | y_{1:T}]}
%      cov_gr:  1xT cell matrix  -  cov(X)_{k|T}  = cov(X_k | y_{1:T})
%                                   (each element is a dxd matrix)
%     cov_wlk:  1xS cell matrix  -  cov(X)_{walk(s)|1:walk(s)}  = cov(X_{walk(s)} | y_{1:walk(s)})
%                                   (each element is a dxd matrix)
% 
% Author: geopapa

    dateStart = datestr(now);
    
    if ~isempty(set) && length(walk)~=length(set)
        error('If a measurement set is specified, it should have the same length with the corresponding visit walk.');
    end
    
    if any(strcmp(varargin,'print')),  printFlag = true;  else  printFlag = false;  end
    
    if ~isempty(set)        
        for k = 1:size(Y,2)
            idx = walk==k;
            if ~isempty(set(idx)) && (max(set(idx)) > length(Y{k}) || min(set(idx)) < 0)
                error('The measurement set for the specified walk cannot have values greater than %d or less than 0 corresponding to the observation set Y{%d}', length(Y{k}), k);
            end
        end
    end
    
    d = size(A{1},1);  % d: dimension of each hidden variable
    T = size(Y,2);     % T: number of time points
    S = length(walk);  % length of visit walk
    
    if ~isempty(expl)        
        for k = 1:T
            if isempty(expl{k})
                expl{k} = 1:length(Y{k});
            end
        end        
    else
        expl = cell(1,T);
        for k = 1:T
            expl{k} = 1:length(Y{k});
        end
    end
    
    % Find up to which point is the walk forward
    fwdIdx = walk(2:end)-walk(1:end-1);
    fwdIdx = find(fwdIdx<0);
    if isempty(fwdIdx)
        fwdIdx = S;
    end
    fwdIdx = fwdIdx(1) -1;
    
    % Store the initial values of mu1, cov1 since they will change later on in the walk loop
    mu1_init  = mu1;
    cov1_init = cov1;
    
    gr_set    = zeros(1,S);
    incr_gr   = zeros(1,S);
    cnd_ig_gr = zeros(1,S);
    
    cov_wlk = [];
    if nargin==16 && flagCovWlk
        cov_wlk = cell(1,S);
    else
        flagCovWlk = false;
    end

    A_gr   = A(1);
    muV_gr = muV(1);
    Q_gr   = Q(1);
    Y_gr   = cell(1);
    C_gr   = cell(1);
    muW_gr = cell(1);
    R_gr   = cell(1);
    c_gr   = 0;
    
    for k = 1:walk(1)-1
        mu1  = A{k}*mu1        + muV{k};
        cov1 = A{k}*cov1*A{k}' +   Q{k};
    end    
    
    updTime_KF = nan(1,S);
    
    for s = 1:S
        if printFlag && (mod(s,10)==0 || s==S || s==1)
            disp(['Iter ', num2str(s), '/', num2str(S), ' +++ Start: ', dateStart, '. Now: ', datestr(now), '.']);
        end        
        
        if ~isempty(set) && set(s) ~= 0
            curr_gr_idx = set(s);            
        else
            curr_gr_idx = findBestEl(s, walk, expl, Y, c, c_gr, budget, lambda, A_gr, muV_gr, Q_gr, C, muW, R, mu1, cov1);
        end

        % Budget is exceeded
        if all(c_gr + c(expl{walk(s)},walk(s)) > budget)
            break;
        end
        
        if curr_gr_idx ~= 0                                    % discard from the exploration set of the current element of the walk
            expl{walk(s)}(expl{walk(s)} == curr_gr_idx) = [];  % the measurement that has just been considered
        else                                                   % or discard all the measurements if no selection has been made
            expl{walk(s)} = [];                                % (since no measurement will be valid for selection later due to submodularity)
        end

        % Update the initial mean, covariance for the next iteration (next item in the walk)
        % (if no item has been greedily selected, the initial mean and cov remain the same)
        if false && curr_gr_idx ~= 0       
            Y_gr{1}   = Y{walk(s)}(curr_gr_idx);
            C_gr{1}   = C{walk(s)}(curr_gr_idx,:);
            muW_gr{1} = muW{walk(s)}(curr_gr_idx);
            R_gr{1}   = R{walk(s)}(curr_gr_idx,curr_gr_idx);
            [~, cov_pr, mean_up, cov_up] = kalm_filt(Y_gr, A_gr, muV_gr, Q_gr, C_gr, muW_gr, R_gr, mu1, cov1, [], [], [], false, false, false);
            
            [max_incr_val, max_cnd_ig]   = eval_incr(cov_pr{1}, cov_up{1}, lambda, c(curr_gr_idx,walk(s)));            
            if max_incr_val < 0,  max_incr_val = 0;  end
            
            gr_set(s)    = curr_gr_idx;
            c_gr         = c_gr + c(curr_gr_idx,walk(s));
            incr_gr(s)   = max_incr_val;
            cnd_ig_gr(s) = max_cnd_ig;
            
            mu1  = mean_up(:,1);   % it would be the updated mean and covariance of the
            cov1 = cov_up{1};      % greedily selected item            
        end
        
        if flagCovWlk              % If an item has been (greedily) selected cov1 is essentially cov_up{1}.
            cov_wlk{s} = cov1;     % Otherwise, it stays the same with previous step
        end                     

        % Next element of the walk would correspond to a different state (time point)
        % Correct initial mean (mu1), covariance (cov1) for that
        if s~=S && walk(s+1)~=walk(s)
            tS = tic;
            if s <= fwdIdx
                for k = walk(s):walk(s+1)-1
                    mu1  = A{k}*mu1        + muV{k};
                    cov1 = A{k}*cov1*A{k}' +   Q{k};
                end
            else
                T_curr = max(walk(1:s+1));
                [Y_gr_cw, C_gr_cw, muW_gr_cw, R_gr_cw] = greedyPrms(T_curr, walk(1:s), gr_set(1:s), Y, C, muW, R);
                [mean_pr, cov_pr, mean_up, cov_up]     = kalm_filt(Y_gr_cw, A(1:T_curr), muV(1:T_curr), Q(1:T_curr), C_gr_cw, muW_gr_cw, R_gr_cw, mu1_init, cov1_init, [], [], [], false, false, true);
                mean_pr                                = mean_pr(:,walk(s+1):end);
                cov_pr                                 = cov_pr(walk(s+1):end);
                mean_up                                = mean_up(:,walk(s+1):end);
                cov_up                                 = cov_up(walk(s+1):end);
                [mean_gr, cov_gr]                      = kalm_smo(A(walk(s+1):T_curr), mean_pr, cov_pr, mean_up, cov_up, false, false);
                mu1                                    = mean_gr(:,1);
                cov1                                   = cov_gr{1};
            end
            tE = toc(tS);
            updTime_KF(s) = tE;
        end
        if s==1 || mod(s,10) == 0,  save(['updTime_KF', num2str(suffx), '.mat'], 's', 'updTime_KF');  end
    end
    
    f_gr  = sum(incr_gr);
    ig_gr = sum(cnd_ig_gr);
end

function curr_gr_idx = findBestEl(s, walk, expl, Y, c, c_gr, budget, lambda, A_gr, muV_gr, Q_gr, C, muW, R, mu1, cov1)    
    Y_gr         = cell(1);
    C_gr         = cell(1);
    muW_gr       = cell(1);    
    R_gr         = cell(1);
    max_incr_val = 0;
    curr_gr_idx  = 0;
    for u = expl{walk(s)}
        if c_gr + c(u,walk(s)) <= budget
            Y_gr{1}   = Y{walk(s)}(u);
            C_gr{1}   = C{walk(s)}(u,:);
            muW_gr{1} = muW{walk(s)}(u);
            R_gr{1}   = R{walk(s)}(u,u);
            
            [~, cov_pr, ~, cov_up] = kalm_filt(Y_gr, A_gr, muV_gr, Q_gr, C_gr, muW_gr, R_gr, mu1, cov1, [], [], [], false, false, false);
            [incr_val, ~]          = eval_incr(cov_pr{1}, cov_up{1}, lambda, c(u,walk(s)));
                        
            if incr_val > max_incr_val
                max_incr_val = incr_val;
                curr_gr_idx  = u;
            end
        end
    end
end

function [Y_gr, C_gr, muW_gr, R_gr] = greedyPrms(T, walk, gr_set, Y, C, muW, R)
    Y_gr   = cell(1,T);
    C_gr   = cell(1,T);
    muW_gr = cell(1,T);
    R_gr   = cell(1,T);
    for k = 1:T
        idx                     = walk==k;
        gr_set_cw               = gr_set(idx);
        gr_set_cw(gr_set_cw==0) = [];
        gr_set_cw               = sort(gr_set_cw);
        Y_gr{k}                 = Y{k}(gr_set_cw);
        C_gr{k}                 = C{k}(gr_set_cw,:);
        muW_gr{k}               = muW{k}(gr_set_cw);
        R_gr{k}                 = R{k}(gr_set_cw,gr_set_cw);
    end
end