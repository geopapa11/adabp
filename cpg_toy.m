function cpg_toy(N_arr,varargin)
% function cpg_toy(N_arr, varargin)
% 
% This function generates an HMM of length N, loads the HMM parameters
% learnt from the CpG Searcher (ran on human chromosome 1).
% It accepts as an argument or produces a mutation sequence and performs a 
% computational mutagenesis analysis by applying the data on RCTreeMP 
% [Sumer et al., 2011], Adaptive MP (AdaMP) or both.
%
% INPUT
%     N_arr:  1xK array    -  Each element represents the length of chain (number of nodes)
%  varargin: cell array    -  Contains pairs of arguments
%                             * If varargin{i} = 'isMeasOrderInc', 
%                               then if varargin{i+1} = true, an increasing  
%                               measurement order will be generated,
%                               if varargin{i+1} = false, a random   
%                               measurement order will be generated.
%                               If flag 'isMeasOrderInc' is not provided, a  
%                               random measurement order will be generated 
%                               by default.
%                             * If varargin{i} = 'onlyadamp', data will be 
%                               run only on AdaMP.
%                             * If varargin{i} = 'onlyrctreemp', data will  
%                               be run only on RCTreeMP [Sumer et al., 2011]
%                               Otherwise, data will be run on both methods.
%                             * If varargin{i} = 'w', then varargin{i+1} 
%                               the variable representing the measurement 
%                               order.
% 
% OUTPUT
% A .mat file with all the important output saved as: 'time_CpG_N_*
% 
% EXAMPLES
%
% One example is if you create an HMM of length 10^5, load a customized
% measurement order from file (that has to be consistent with the length 
% of HMM) and run it only on the RCTreeMP algorithm by [Sumer et al., 2011]
% E.g.,
% N = 10^5;
% w = [3 23 21 3 1001 22001 57];
% cpg_toy(N, 'w', w, 'onlyrctreemp');
% 
% Another example is generating an HMM of length 10^4 with an increasing 
% measurement order which runs on both methods (RCTreeMP, AdaMP).
% Namely:
% N = 10^4;
% cpg_toy(N,'isMeasOrderInc',true);
% 
% If you wish to reproduce Fig. 4a of the paper, run the following command:
% N_arr = 10.^[2 3 4 5];
% cpg_toy(N_arr,'isMeasOrderInc',true);
% (
%   Unfortunately, we noticed that RCTreeMP might crash in newer versions 
%   of MATLAB. If you wish to run AdaMP alone, please type the command:
%   N_arr = 10.^[2 3 4 5];
%   cpg_toy(N_arr,'isMeasOrderInc',true,'onlyadamp');
% )
% 
% If you wish to reproduce Fig. 4b of the paper, run the following command:
% N_arr = 10.^[2 3 4 5];
% cpg_toy(N_arr);
% (
%   Unfortunately, we noticed that RCTreeMP might crash in newer versions 
%   of MATLAB. If you wish to run AdaMP alone, please type the command:
%   N_arr = 10.^[2 3 4 5];
%   cpg_toy(N_arr,'onlyadamp');
% )
% 
% Author: geopapa
% $ Date: 2015/07/01 22:43:11 $

close all;

% CpGIslandSearcher.mat contains the following variables:
%     -     pi: initialization matrix ( pi(1): Propability of the first nucleotide to not correspond to a CpG island )
%     -      a: transition matrix ( a(i,j): Propability of transitioning from i to j )
%                              i = 1 => 0, i = 2 => 1 (CpG position)
%     -      e: emission matrix ( 
%                            e(1,1): Propability of emitting A if not a CpG
%                            e(1,2): Propability of emitting C if not a CpG
%                            e(1,3): Propability of emitting G if not a CpG
%                            e(1,4): Propability of emitting T if not a CpG
%                            e(2,1): Propability of emitting A if a CpG
%                            e(2,2): Propability of emitting C if a CpG
%                            e(2,3): Propability of emitting G if a CpG
%                            e(2,4): Propability of emitting T if a CpG
%   -      CpG: Prediction by CpG Island Searcher
%   -      seq: Nucleotide sequence (as sequence of numbers)
%   -   seqstr: Nucleotide sequence (as string of A, C, G, T characters)
%   - startCpG: Start of CpG island
%   -   endCpG: End of CpG island
load('CpGIslandSearcher.mat');

clear seqstr CpG;

disp_intv = 1;
[card, obs_size] = size(e);
K = length(N_arr);

time1         = nan(size(N_arr));
time2         = nan(size(N_arr));
time1_upd_inn = cell(size(N_arr));
time2_upd_inn = cell(size(N_arr));
time1_map_inn = cell(size(N_arr));
time2_map_inn = cell(size(N_arr));
path_len_ww   = cell(size(N_arr));
avg_path_len  = zeros(size(N_arr));
mapseq1_diff  = cell(size(N_arr));
mapseq2_diff  = cell(size(N_arr));

if any(strcmp(varargin,'isMeasOrderInc'))
    isMeasOrderInc = varargin{find(strcmp(varargin,'isMeasOrderInc'))+1};
else
    isMeasOrderInc = false;
end

if any(strcmp(varargin,'w'))
    w = varargin{find(strcmp(varargin,'w'))+1};
    if all(w(2:end)-w(1:end-1) >= 0), isMeasOrderInc = true;  else  isMeasOrderInc = false;  end
else
    w = [];
end

if any(strcmp(varargin,'onlyrctreemp'))
    runrctreemp = true;
    runadamp    = false;
elseif any(strcmp(varargin,'onlyadamp'))
    runrctreemp = false;
    runadamp    = true;
else
    runrctreemp = true;
    runadamp    = true;
end

if isMeasOrderInc,  sfx = '_inc';  else  sfx = '';  end

disp(['Is measurement order w increasing?: ', num2str(isMeasOrderInc)]);

if runrctreemp && runadamp
    outfile = ['time_CpG_N_',num2str(N_arr(end)), sfx, '.mat'];
elseif runrctreemp
    outfile = ['time_CpG_N_',num2str(N_arr(end)), sfx, '_rctreemp.mat'];
else
    outfile = ['time_CpG_N_',num2str(N_arr(end)), sfx, '_adamp.mat'];
end

%% Observed sequence
y_all = seq;

%% Construct the potentials
N = N_arr(end);

% Node potentials
phi_all           = cell(1,N);
phi_all(y_all>0)  = num2cell(e(:,y_all(y_all>0)),1);
phi_all(y_all==0) = {ones(card,1)};
phi_all{1}        = pi.*phi_all{1};

clear y_all seq;

% Pairwise potential (common across all pairs)
psi = a;

%% Construct the Markov chain
E_all = [(1:N-1); 2:N];
E_all = E_all';

if numel(N_arr) == 1
    N       = N_arr(1);
    phi_all = phi_all(1:N);
    E_all   = E_all(1:N-1,:);
end

dateStart = datestr(now);
t_s = tic;
r = 1;
n = 1;
for N = N_arr
    disp(['N(', num2str(n), ')=', num2str(N)]);
    
    Ne = N-1;
    
    %% Isolate a part of the whole sequence
    phi  = phi_all(1:N);
    E    = E_all(1:N-1,:);
    root = randi(1);
    
    %% Construct mutation sequence, w if it is not provided as an argument
    if isempty(w)
        if N <= 10000
            if isMeasOrderInc,  w = 1:N;    else  w = randperm(N);                       end
        else
            if isMeasOrderInc,  w = 1:2:N;  else  w = randperm(N); w = w(1:round(N/2));  end
        end
    end
    
    M = length(w);
        
    %% Choose mutations
    y = randi(obs_size,1,N);
    
    %% Compare timing for RCTreeMP [Sumer et al., 2011]
    time1_upd_inn{n} = zeros(1,M-1);
    time1_map_inn{n} = zeros(1,M-1);
    mapseq1_diff{n}  = zeros(1,M-1);
	phi1             = phi;
    A                = sparse([1:N, E(:)'], [1:N, repmat(N+1:N+Ne,1,2)], 1);  % Adjacent matrix
    psi_tmp          = cell(1,Ne);
    psi_tmp(:)       = {a};
    F                = [phi, psi_tmp];
    t_s1       = tic;
    time1_init = toc(t_s1);
    time1(r,n) = toc(t_s1);
    
    if runrctreemp
        % Build RC-tree structure for Max-Product
        T = rctreeMP(A,F,8);
        time1_init = toc(t_s1);
        j = 1;
        if w(j) == 1,  phi1{w(j)} = pi.*e(:,y(j));  else  phi1{w(j)} = e(:,y(j));  end
        T            = updateFactorMP(T, w(j), phi1{w(j)});
        mapseq1_prv  = double(T.map);
        for j = 2:M
            t_s1_upd = tic;
            if w(j) == 1,  phi1{w(j)} = pi.*e(:,y(j));  else  phi1{w(j)} = e(:,y(j));  end
            [T, t_map] = updateFactorMP(T, w(j), phi1{w(j)});
            
            time1_upd_inn{n}(j-1) = toc(t_s1_upd);
            time1_upd_inn{n}(j-1) = time1_upd_inn{n}(j-1) - t_map;
            
            t_s1_map = tic;
            mapseq1 = double(T.map);
            time1_map_inn{n}(j-1) = toc(t_s1_map);
            time1_map_inn{n}(j-1) = time1_map_inn{n}(j-1) + t_map;
            
            mapseq1_diff{n}(j-1) = sum(mapseq1~=mapseq1_prv);
            mapseq1_prv          = mapseq1;
            
            if mod(j,10000) == 0
                t_e = toc(t_s);
                disp(['    RCTreeMP :: Update # ', num2str(j), '/', num2str(M), ' has been completed. Start: ', dateStart, '. Elapsed time (', datestr(now), '): ', num2str(t_e), ' sec.']);
                save(['time_CpG_N_',num2str(N_arr(1)), sfx, '_rctree.mat'],'w','card','N_arr','time1','time1_init','time1_upd_inn','time1_map_inn','mapseq1_diff');
            end
        end
        time1(r,n) = toc(t_s1);
    end
    
    clear phi1;
    
    %% Compare timing for Adaptive MP; AdaMP
    time2_upd_inn{n} = zeros(1,M-1);
    time2_map_inn{n} = zeros(1,M-1);
    mapseq2_diff{n}  = zeros(1,M-1);
    path_len_ww{n}   = zeros(1,M-1);
    phi2             = phi;
    t_s2       = tic;
    time2_init = toc(t_s2);
    time2(r,n) = toc(t_s2);
    
    if runadamp
        % Build AdaMP structure for Max-Product
        abp = AdaBP(E,phi,psi,'root',root,'max');
        time2_init = toc(t_s2);
        j = 1;
        if w(j) == 1,  phi2{w(j)} = pi.*e(:,y(j));  else  phi2{w(j)} = e(:,y(j));  end
        abp.setNodePot(w(j), phi2{w(j)});
        mapseq2_prv = abp.mapseq();
        for j = 2:M
            t_s2_upd = tic;
            if w(j) == 1,  phi2{w(j)} = pi.*e(:,y(j));  else  phi2{w(j)} = e(:,y(j));  end
            abp.setNodePot(w(j), phi2{w(j)});
            plw = abp.propagate(w(j-1), w(j), true);
            time2_upd_inn{n}(j-1) = toc(t_s2_upd);
            
            t_s2_map = tic;
            mapseq2 = abp.mapseq();
            time2_map_inn{n}(j-1) = toc(t_s2_map);
            
            path_len_ww{n}(j-1) = plw;
            
            mapseq2_diff{n}(j-1) = sum(mapseq2~=mapseq2_prv);
            mapseq2_prv          = mapseq2;
            
            if mod(j,10000) == 0
                t_e = toc(t_s);
                disp(['    AdaMP :: Update # ', num2str(j), '/', num2str(M), ' has been completed. Start: ', dateStart, '. Elapsed time (', datestr(now), '): ', num2str(t_e), ' sec.']);
                save(['time_CpG_N_',num2str(N_arr(1)), sfx, '_adabp.mat'],'w','card','N_arr','time2','time2_init', 'time2_upd_inn','time2_map_inn', 'path_len_ww','mapseq2_diff');
            end
        end
        time2(r,n) = toc(t_s2);
    end
    
    avg_path_len(r,n) = mean(path_len_ww{n});
    
    if mod(r,disp_intv) == 0
        t_e = toc(t_s);
        disp(['    Sample # ', num2str(r), '/1 has been processed. Start: ', dateStart, '. Elapsed time (', datestr(now), '): ', num2str(t_e), ' sec.']);
    end
    n = n + 1;
    save(outfile,'w','card','N_arr','time1','time2','time1_init','time1_upd_inn','time1_map_inn', 'time2_init', 'time2_upd_inn','time2_map_inn', 'path_len_ww','avg_path_len','mapseq1_diff','mapseq2_diff');
    
    if ~any(strcmp(varargin,'w')),  w = [];  end
end

rt_time1_to_time2_emp = time1./time2;

save(outfile,'rt_time1_to_time2_emp','-append');

disp('');
disp(['card = [',  num2str(card),  ']']);
disp(['N = [',     num2str(N_arr), ']']);
disp('');
disp(['Empirical avg path len = [', num2str(avg_path_len), ']']);
disp('');
disp(['Empirical Ratio Time(RCTreeMP)/Time(AdaMP) = [', num2str(rt_time1_to_time2_emp), ']']);
disp('');

h1 = figure;
[hAx,hBar,hLine] = plotyy(log10(N_arr), time1./time2, [log10(N_arr)', log10(N_arr)'], [time1', time2'], 'bar','semilogy');
set(hBar, 'FaceColor', 'k', 'BarWidth', .4);
set(hAx,{'ycolor'},{'k';'b'});
set(hLine(1), 'Color', 'b', 'LineWidth', 2);
set(hLine(2), 'Color', 'g', 'LineWidth', 2);
set(hAx,'FontSize',30);
xlabel('$N$', 'interpreter','latex');
legend(hLine, {'RCTreeMP', 'AdaMP'},'FontSize',26, 'Location','NorthWest');
set(hAx, 'XTickLabel', []);
xticks = get(hAx, 'XTick');
xticks = xticks{1};
arrayfun(@(x)text(x-.1,-.55,texlabel(sprintf('10^%d',x)), 'FontSize', 30),xticks,'UniformOutput',false);
axis tight;
saveas(h1,'cpg_toy_time.fig');

h2 = figure;
for i = 1:K
    path_len  = path_len_ww{i};
    time1_inn = time1_upd_inn{i} + time1_map_inn{i};
    time2_inn = time2_upd_inn{i} + time2_map_inn{i};
    subplot(K, 1, i);
    semilogx(path_len, time1_inn./time2_inn, '.', 'Color', 'k', 'Markersize', 14);
    hold on;
    semilogx(1:N_arr(end), ones(1,N_arr(end)), 'r');
    if i~=K, set(gca, 'XTickLabel', []);  end
    set(gca,'XLim', [1 N_arr(end)]);
    set(gca,'YScale', 'log');
    set(gca,'FontSize',30);
    if i == 1,  title('Speedup of AdaMP over RCTreeMP', 'FontSize', 30);            end
    if i == K,  xlabel('Distance between consecutive w elements', 'FontSize', 30);  end
end
saveas(h2,'cpg_toy_hist.fig');