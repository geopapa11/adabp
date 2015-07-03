function [gr_set, f_gr, ig_gr, incr_gr, cnd_ig_gr, mean_gr, log_cov_gr, log_cov_wlk] = greedy_time_upd_adabp(suffx, walk, set, expl, Y, c, budget, lambda, A, muV, Q, C, muW, R, mu1, cov1, varargin)
% function [gr_set, f_gr, ig_gr, incr_gr, cnd_ig_gr, mean_gr, log_cov_gr, log_cov_wlk] = greedy_time_upd_adabp(suffx, walk, set, expl, Y, c, budget, lambda, A, muV, Q, C, muW, R, mu1, cov1, varargin)
%
% This function finds the greedy set, its reward, MI, the incremental 
% reward and conditional MI at each step for the specified visit walk.
%
% Lastly, it gives the greedy estimates of the parameter of interest:
% E(X)_{k|T}   = { E[X_1 | y_{1:T}], E[X_2 | y_{1:T}], ..., E[X_T | y_{1:T}] }
% and
% log(det(cov(X)_{k|T})) = { log(det(cov(X_1 | y_{1:T}))), log(det(cov(X_2 | y_{1:T}))), ..., log(det(cov(X_T | y_{1:T}))) }.
%
% This function uses Adaptive BP to find the covariance of the next element 
% of the walk (it does not use Kalman filtering or smoothing).
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
%  varargin: cell array            -  Contains optional arguments
%                                     * If varargin{i} = 'flagCovWlk',
%                                       then we store the cov updates for each step of the walk
%                                       (Default value: Cov values are not stored)
%
%
% OUTPUT
%
%       gr_set:  1xS  array      -  Greedy set showing which measurements have been
%                                   chosen greedily for each point in the walk
%         f_gr:  1x1  scalar     -  Reward of the greedy solution
%        ig_gr:  1x1  scalar     -  MI of the greedy solution
%      incr_gr:  1xT  array      -  Incremental reward of the greedy solution for each time point
%                                   incr_gr(k) = I(X; Y_k | y_{1:k-1}) - lambda*c_k
%    cnd_ig_gr:  1xT  array      -  Conditional MI of the greedy solution for each time point
%                                   cnd_ig_gr(k) = I(X; Y_k | y_{1:k-1})
%      mean_up:  1xL cell array  -  E(X)_{k|k}    = {E[x_1 | y_{1}], E[x_2 | y_{1:2}], ..., E[x_T | y_{1:T}]}
%                                   Each element is a dxT  matrix
%                                   representing the updated means for fixed
%                                   values of A, Q, R
%   log_cov_up:  LxT array       -  log(det(cov(X)_{k|k}))  = log(det(cov(X_k | y_{1:k})))
%                                   (The element of each row is the log
%                                   determinant of the covariance for fixed
%                                   values of A, Q, R
%     mean_gr:  1xL cell array   -  E(X)_{k|T}    = {E[x_1 | y_{1:T}], E[x_2 | y_{1:T}], ..., E[x_T | y_{1:T}]}
%                                   Each element is a dxT  matrix
%                                   representing the smoothed means for fixed
%                                   values of A, Q, R
%  log_cov_up:  LxT array        -  log(det(cov(X)_{k|T}))  = log(det(cov(X_k | y_{1:T})))
%                                   (The element of each row is the log
%                                   determinant of the covariance for fixed
%                                   values of A, Q, R
% log_cov_wlk:  1xS array        -  log(det(cov(X)_{walk(s)|1:walk(s)})) = log(det(cov(X_{walk(s)} | y_{1:walk(s)})))
%
% Author: geopapa
% $ Date: 2014/01/22 21:03:42 $

    dateStart = datestr(now);
    
    if ~isempty(set) && length(walk)~=length(set)
        error('If a measurement set is specified, it should have the same length with the corresponding visit walk.');
    end

    if ~isempty(set)
        for k = 1:size(Y,2)
            idx = walk==k;
            if ~isempty(set(idx)) && (max(set(idx)) > length(Y{k}) || min(set(idx)) < 0)
                error('The measurement set for the specified walk cannot have values greater than %d or less than 0 corresponding to the observation set Y{%d}', length(Y{k}), k);
            end
        end
    end

    gpuExists = gpuDeviceCount > 0;

    if ~any(strcmp(varargin,'useGPU'))
        useGPU = true;  % default is to use GPU if it exists
    else
        useGPU = varargin{find(strcmp(varargin,'useGPU'))+1};
    end

    T = size(Y,2);     % number of time points
    S = length(walk);  % length of visit walk

    % Set the exploration sets wherever they are empty
    expl = set_expl(expl, Y);

    % Check if we want to store the covariance after the greedy
    % incorporation of a measurement and set some parameters
    [flagCovWlk, log_cov_wlk, printFlag, gr_set, incr_gr, cnd_ig_gr, c_gr] = ...
                                  set_prms(S, gpuExists, useGPU, varargin);

    % Store the initial values of mu1, cov1 since they will change later on in the walk loop
    mu1_init = mu1;  cov1_init = cov1;

    % Determine the variables that each measurement is drawn from
    Ic = cellfun(@(x) ne(x,0), C, 'UniformOutput', false);
    
    % Initialize Adaptive BP
    [abp, mu1, cov1] = init_adabp(walk, A, muV, Q, Y, mu1_init, cov1_init, gpuExists, useGPU, T);

    updTime_ABP = nan(1,S);
    
    for s = 1:S
        if printFlag && (mod(s,10)==0 || s==S || s==1)
            disp(['Iter ', num2str(s), '/', num2str(S), ' +++ Start: ', dateStart, '. Now: ', datestr(now), '.']);
        end

        expl_cur = expl{walk(s)};  Y_cur = Y{walk(s)};    c_cur = c(:,walk(s));
        C_cur    = C{walk(s)};  muW_cur  = muW{walk(s)};  R_cur = R(:,walk(s));
        Ic_cur   = Ic{walk(s)};

        if ~isempty(set) && set(s) ~= 0
            gr_idx = set(s);
        else
            [max_incr, gr_idx] = findBestEl(expl_cur, Y_cur, c_cur, c_gr, budget, lambda, C_cur, muW_cur, R_cur, mu1, cov1, Ic_cur, gpuExists, useGPU);
        end

        % Budget is exceeded
        if all(c_gr + c_cur(expl{walk(s)}) > budget)
            break;
        end
        
        % Store the results obtained in the current greedy step
        if false && gr_idx ~= 0
            gr_set(s)  = gr_idx;    c_gr         = c_gr + c_cur(gr_idx);
            incr_gr(s) = max_incr;  cnd_ig_gr(s) = max_incr + lambda*c_cur(gr_idx);

            % Remove from the current exploration set the measurement that has just been obtained
            expl{walk(s)}(expl{walk(s)} == gr_idx) = [];
        else
            % Due to submodularity, all the remaining measurements are guaranteed
            % to be of no value, if no selection has been made in this step
            expl{walk(s)} = [];
        end

        % Find mu1, cov1 of next walk element
        if s~=S  && gr_idx ~= 0
            tS = tic;
            wlk_cur = walk(s);
            wlk_nxt = walk(s+1);
            Y_u     = Y_cur(gr_idx);
            C_u     = C_cur(gr_idx,:);
            muW_u   = muW_cur(gr_idx);
            R_u     = cell(size(R_cur));
            for r = 1:length(R_cur)
                R_u{r} = R_cur{r}(gr_idx,gr_idx);
            end
            [mu1, cov1, log_cov_wlk(s)] = upd_nxt_mu_cov(abp, wlk_cur, wlk_nxt, Y_u, C_u, muW_u, R_u, gpuExists, useGPU, flagCovWlk);
            tE = toc(tS);
            updTime_ABP(s) = tE;
        end
        
        if s==1 || mod(s,10) == 0,  save(['updTime_AdaBP', num2str(suffx), '.mat'], 's', 'updTime_ABP');  end
    end  % end of for s = 1:S

    % Find the means and log(det(cov)) of all the hidden states after the end of the greedy algorithm
    [mean_gr, log_cov_gr] = eval_mu_cov(abp, T, gpuExists, useGPU);
    
    f_gr  = sum(incr_gr);
    ig_gr = sum(cnd_ig_gr);
end  % end of greedy method

function expl = set_expl(expl, Y)
    T = size(Y,2);
    if ~isempty(expl)
        for k = 1:T,  if isempty(expl{k}),  expl{k} = 1:length(Y{k});  end
        end
    else
        expl = cell(1,T);
        for k = 1:T,  expl{k} = 1:length(Y{k});  end
    end
end

function [flagCovWlk, log_cov_wlk, printFlag, gr_set, incr_gr, cnd_ig_gr, c_gr] = ...
        set_prms(S, gpuExists, useGPU, varargin)
    % Check if output will be printed (default is no)
    if any(strcmp(varargin{:},'print')),  printFlag = true;  else  printFlag = false;  end

    % Check if we want to store the covariance after the greedy incorporation of a measurement
    if gpuExists && useGPU
        if any(strcmp(varargin{:},'flagCovWlk'))
            flagCovWlk = gpuArray.true;  log_cov_wlk = gpuArray.zeros(1,S,'single');
        else
            flagCovWlk = gpuArray.false; log_cov_wlk = [];
        end
        gr_set    = gpuArray.zeros(1,S,'single');
        incr_gr   = gpuArray.zeros(1,S,'single');
        cnd_ig_gr = gpuArray.zeros(1,S,'single');
        c_gr      = gpuArray.zeros(1,1,'single');
    else
        if any(strcmp(varargin{:},'flagCovWlk'))
            flagCovWlk = true;  log_cov_wlk = zeros(1,S);
        else
            flagCovWlk = false; log_cov_wlk = [];
        end
        gr_set    = zeros(1,S,'single');
        incr_gr   = zeros(1,S);
        cnd_ig_gr = zeros(1,S);
        c_gr      = 0;
    end
end

function [abp, mu1, cov1] = init_adabp(walk, A, muV, Q, Y, mu1, cov1, gpuExists, useGPU, T)
    if iscell(A),  L = size(A,1);  else  L = 1;  end    
    abp = cell(1,L);
    
    if ~iscell(mu1),   tmp = cell(1,L);  tmp(:) = {mu1};   mu1  = tmp;  end
    if ~iscell(cov1),  tmp = cell(1,L);  tmp(:) = {cov1};  cov1 = tmp;  end
    clear tmp;
    
    wlk_first_el = walk(1);
    if gpuExists && useGPU
        for r = 1:L
            if iscell(A),                         A_r   = A(r,:);    else  A_r    = A;     end
            if length(A_r)~=T && length(A_r)==1,  A_r   = A_r{1};    end
            if iscell(muV)  && size(muV,1)==L,    muV_r = muV(r,:);  else  muV_r  = muV;   end
            if iscell(Q)    && size(Q,1)==L,      Q_r   = Q(r,:);    else  Q_r    = Q;     end
            if length(Q_r)~=T && length(Q_r)==1,  Q_r   = Q_r{1};    end            
            abp{r}  = AdaBP_HMM(A_r, muV_r, Q_r, Y, mu1{r}, cov1{r});
            [h, J]  = abp{r}.eval_mrg(wlk_first_el);
            cov1{r} = J\eye(size(J));
            cov1{r} = (cov1{r} + cov1{r}')/2;
            mu1{r}  = cov1{r}*h;
        end
    else
        for r = 1:L
            if iscell(A),                         A_r   = A(r,:);    else  A_r    = A;     end
            if length(A_r)~=T && length(A_r)==1,  A_r   = A_r{1};    end
            if iscell(muV)  && size(muV,1)==L,    muV_r = muV(r,:);  else  muV_r  = muV;   end
            if iscell(Q)    && size(Q,1)==L,      Q_r   = Q(r,:);    else  Q_r    = Q;     end
            if length(Q_r)~=T && length(Q_r)==1,  Q_r   = Q_r{1};    end
            abp{r}  = AdaBP_HMM(A_r, muV_r, Q_r, Y, mu1{r}, cov1{r});
            [h, J]  = abp{r}.eval_mrg(wlk_first_el);
            cov1{r} = J\eye(size(J));
            cov1{r} = (cov1{r} + cov1{r}')/2;
            mu1{r}  = cov1{r}*h;
        end
    end
end

function [mu1, cov1, log_cov_wlk_cur] = upd_nxt_mu_cov(abp, wlk_cur, wlk_nxt, Y_cur, C_cur, muW_cur, R_cur, gpuExists, useGPU, flagCovWlk)
    mu1  = cell(size(abp));
    cov1 = cell(size(abp));
    log_cov_wlk_cur = NaN;

    if gpuExists && useGPU
        log_cov_wlk_arr = gpuArray.zeros(size(abp),'single');
    else
        log_cov_wlk_arr = zeros(size(abp));
    end
    
    for r = 1:length(abp)
        abp{r}.update(wlk_cur, Y_cur, C_cur, muW_cur, R_cur{r});
    end
    
    if flagCovWlk
        for r = 1:length(abp)
            [~, J]   = eval_mrg(obj, wlk_cur);
            cov1_tmp = J\eye(size(J));
            cov1_tmp = (cov1_tmp + cov1_tmp')/2;
            log_cov_wlk_arr(r) = sum(log(eig(cov1_tmp)));
        end
        log_cov_wlk_cur = nanmean(log_cov_wlk_arr);
    end
    
    if gpuExists && useGPU
        for r = 1:length(abp)
            abp{r}.propagate(wlk_cur, wlk_nxt);
            [h, J]  = abp{r}.eval_mrg(wlk_nxt);
            cov1{r} = J\eye(size(J));
            cov1{r} = (cov1{r} + cov1{r}')/2;
            mu1{r}  = cov1{r}*h;
        end
    else
        if length(abp) > 1
            parfor r = 1:length(abp)
                abp{r}.propagate(wlk_cur, wlk_nxt);
                [h, J]  = abp{r}.eval_mrg(wlk_nxt);
                cov1{r} = J\eye(size(J));
                cov1{r} = (cov1{r} + cov1{r}')/2;
                mu1{r}  = cov1{r}*h;
            end
        else
            for r = 1:length(abp)
                abp{r}.propagate(wlk_cur, wlk_nxt);
                [h, J]  = abp{r}.eval_mrg(wlk_nxt);
                cov1{r} = J\eye(size(J));
                cov1{r} = (cov1{r} + cov1{r}')/2;
                mu1{r}  = cov1{r}*h;
            end
        end
    end
end

function [mean_gr, log_cov_gr] = eval_mu_cov(abp, T, gpuExists, useGPU)
% we assume that all hidden variables have the same dimension
    L       = length(abp);
    mean_gr = cell(1,L);
    for r = 1:L
        h          = abp{r}.eval_mrg(1);
        mean_gr{r} = zeros(length(h),T);
    end
    
    log_cov_gr = zeros(L,T);
    
    if gpuExists && useGPU
        for r = 1:L
            abp{r}.reset_msg();
            abp{r}.propagate(1,T);
            abp{r}.propagate(T,1);
        end
    else
        if L > 1
            parfor r = 1:L
                abp{r}.reset_msg();
                abp{r}.propagate(1,T);
                abp{r}.propagate(T,1);
            end
        else
            for r = 1:L
                abp{r}.reset_msg();
                abp{r}.propagate(1,T);
                abp{r}.propagate(T,1);
            end
        end
    end
    
    for r = 1:L
        for k = 1:T
            [h, J]          = abp{r}.eval_mrg(k);
            mean_gr{r}(:,k) = J\h;
            log_cov_gr(r,k) = -sum(log(eig(J)));
        end
    end
end