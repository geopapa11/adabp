function berk_toy(flag, berkData_file, sensCoords_file, bad_sens_idx_file, learnedParams_file, range, theta)
% function berk_toy(flag, berkData_file, sensCoords_file, bad_sens_idx_file, learnedParams_file, range, theta)
% 
% INPUT
%              flag:  scalar      -  Flag denoting which method will be used
%                                    (flag = 0: Kalman filter, flag = 1: Gaussian AdaBP)
%     berkData_file:  cell array  -  Contains pairs of arguments
%   sensCoords_file:  string      -  Sensor locations file
% bad_sens_idx_file:  string      -  List of all sensors with erroneously reported measurements
% learnedParams_file:  mat file   -  Learned parameters by using the sensor measurements
%             range:  array       -  Range of time points under consideration
%             theta:  scalar      -  Bias of the mean of matrix A
% 
% OUTPUT
% A .mat file with all the important output saved as: 'updTime_AdaBP*' run 
% on AdaBP or a .mat file saved as: 'updTime_KF*' run on Kalman
% filter/smoother.
% 
% EXAMPLES
%
% An example is
% flag               = 1; % flag = 1: greedy (with AdaBP), flag = 0: greedy (with Kalman filter)
% berkData_file      = 'data/berkData.mat';
% sensCoords_file    = 'data/sens_coords.txt';
% bad_sens_idx_file  = 'data/bad_sens_idx.txt';
% learnedParams_file = 'data/learnedParams.mat';
% range              = 1:360;
% theta              = 0.4;
% berk_toy(flag, berkData_file, sensCoords_file, bad_sens_idx_file, learnedParams_file, range, theta);
% 
% flag = 0; berkData_file = 'data/berkData.mat'; sensCoords_file = 'data/sens_coords.txt'; bad_sens_idx_file = 'data/bad_sens_idx.txt'; learnedParams_file = 'data/learnedParams.mat'; range = 1:360; theta = .4; berk_toy(flag, berkData_file, sensCoords_file, bad_sens_idx_file, learnedParams_file, range, theta);
% flag = 1; berkData_file = 'data/berkData.mat'; sensCoords_file = 'data/sens_coords.txt'; bad_sens_idx_file = 'data/bad_sens_idx.txt'; learnedParams_file = 'data/learnedParams.mat'; range = 1:360; theta = .4; berk_toy(flag, berkData_file, sensCoords_file, bad_sens_idx_file, learnedParams_file, range, theta);
% 
% Author: geopapa
% $ Date: 2015/07/02 20:12:23 $

    % Import the sensor locations (coordinates)
    formatSpec = '%f%f%f%[^\n\r]';
    delim      = '\t';
    startRow   = 2;
    fid        = fopen(sensCoords_file,'r');
    dataArray  = textscan(fid, formatSpec, 'Delimiter', delim, 'HeaderLines',startRow-1, 'ReturnOnError', false);    
    fclose(fid);
    sens_coords = [dataArray{:, 2} dataArray{:, 3}];
    
    % Import the valid time points file
    badsensidx = [];
    if ~isempty(bad_sens_idx_file)
        formatSpec = '%f%[^\n\r]';
        delim      = '';
        startRow   = 2;
        fid        = fopen(bad_sens_idx_file,'r');
        dataArray  = textscan(fid, formatSpec, 'Delimiter', delim, 'HeaderLines',startRow-1, 'ReturnOnError', false);
        fclose(fid);
        badsensidx = dataArray{:, 1};
    end
    
    load(berkData_file);
    
    mu1_knownpart = Y{5432};
    
    Y = Y(range);
    
    % Convert the observations in double matrix
    ys = cell2mat(Y);
    ys(:,1:58) = [];
    
    % Remove the bad sensors from the sensor list
    sens_coords(badsensidx,:)   = [];
    mu1_knownpart(badsensidx,:) = [];
    ys(badsensidx,:)            = [];
    
    grid_x = unique(sens_coords(:,1))';
    grid_y = unique(sens_coords(:,2))';    
    
    [m, T] = size(ys);       % m: measurement dimension, T: # time points
    d      = length(grid_x)*length(grid_y);  % dimension of X_k
    Y      = num2cell(ys,1);
    
    % Get the map from locations to variable indices
    keys = cell(d,1);
    vals = zeros(d,1);
    cnt  = 1;
    for j = grid_y
        for i = grid_x
           keys{cnt} = ['(', num2str(i), ',', num2str(j), ')'];
           vals(cnt) = cnt;
           cnt       = cnt + 1;
        end
    end
    sub2idx = containers.Map(keys,vals);
    
    % Get the map from variable indices to locations
    idx2sub = zeros(d,2);
    cnt     = 1;
    for j = grid_y
        for i = grid_x
           idx2sub(cnt,:) = [i  j];
           cnt            = cnt + 1;
        end
    end    
    
    % Initialize the parameters
    % Construct M (mean of A)
    M = zeros(d);
    for i = 1:d-1
        for j = i+1:d
            dist   = norm(idx2sub(i,:) - idx2sub(j,:));
            M(i,j) = exp(-theta*dist);
        end
    end
    M = eye(d) + M + M';
    M(abs(M)<.04) = 0;
    
    % Construct C
    C_fxd = zeros(m,d);
    for k = 1:m
        i = sens_coords(m,1);
        j = sens_coords(m,2);
        C_fxd(k,sub2idx(['(', num2str(i), ',', num2str(j), ')'])) = 1;
    end
    
    C    = cell(1,T);
    C(:) = {C_fxd};
    Inan = isnan(ys);
    for t = 1:T
        C{t}(Inan(:,t),:) = zeros(sum(Inan(:,t)),d);
    end
        
    % Construct mu1
    mu1 = zeros(d,1);
    idx = zeros(m,1);
    for s = 1:m
        i      = sens_coords(s,1);
        j      = sens_coords(s,2);
        idx(s) = sub2idx(['(', num2str(i), ',', num2str(j), ')']);
    end    
    mu1(idx) = mu1_knownpart;
    mu1(isnan(mu1)) = 0;
    mu1  = M*mu1;
    cov1 = .1*eye(d);
        
    c      = zeros(m,T);
    budget = Inf;
    lambda = 0;
    
    load(learnedParams_file);
    d   = size(A{1},1);
	muV = zeros(d,1);
	muW = zeros(m,1);
	muV = {muV};
	muW = {muW};
	muV = repmat(muV,1,T);
	muW = repmat(muW,1,T);
    clearvars -except flag Inan T Y c budget lambda A Q C R mu1 cov1 muV muW;
    load('data/walk_arr.mat');  % load predefined walks (each row is a different walk)
    % walk = randperm(T);  % you can define any random walk
    
    for q = 1:size(walk_arr,1)
        disp(['Walk # ', num2str(q), '/', num2str(size(walk_arr,1)), ' is being considered +++ Start: ', datestr(now), '.']);
        
        walk = walk_arr(q,:);

        expl = cell(1,T);
        for t = 1:T
            expl{t} = find(~Inan(:,t));
            if ~isrow(expl{t}),  expl{t} = expl{t}';  end
        end    
        idx = ones(1,T);
        set = zeros(size(walk));
        for s = 1:size(walk,2)
            ii     = randperm(length(expl{walk(s)}),1);
            set(s) = expl{walk(s)}(idx(walk(s)));        
            idx(walk(s)) = idx(walk(s)) + 1;
        end
        
        if flag==1
            [gr_set, f_gr, ig_gr, incr_gr, cnd_ig_gr, mean_gr]             = greedy_time_upd_adabp(q, walk, set, expl, Y, c, budget, lambda, A, muV, Q, C, muW, R, mu1, cov1, 'useGPU', false, 'print');
        else            
            [gr_set2, f_gr2, ig_gr2, incr_gr2, cnd_ig_gr2, ~, ~, mean_gr2] = greedy_time_upd_kf(q, walk, set, expl, Y, c, budget, lambda, repmat(A,1,T), muV, repmat(Q,1,T), C, muW, R, mu1, cov1, [], 'print');
        end
    end
end