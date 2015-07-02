function adabp_toy(card,varargin)
% function adabp_toy(card,varargin)
% 
% This function is a toy example of RCTreeBP [Sumer et al., 2011] 
% and AdaBP performance on synthetically generated trees.
%
% INPUT
%        card:  1x1   scalar  -  Common cardinality among all nodes
%  varargin: cell array       -  Contains pairs of arguments
%                                * If varargin{i} = 'speedup_rt', then varargin{i+1}
%                                  is a 1x1 scalar representing the
%                                  enforced speedup ratio of AdaBP over RCTreeBP
%                                  (this number would affect the composition
%                                  of the measurement order w).
%                                  (default: NaN - no speedup ratio enforced)
%                                * If varargin{i} = 'unbal_rt', then varargin{i+1}
%                                  is a 1x1 scalar representing the ratio
%                                  between 0 and 1 of how unbalanced we
%                                  want the tree to be.
%                                  The higher unbal_rt is, the more
%                                  unbalanced the tree will be (larger
%                                  diameter)
%                                  (default: 0.7)
%                                * If varargin{i} = 'scenario', then varargin{i+1}
%                                  is a string denoting the scenario that
%                                  will be considered.
%                                  If varargin{i+1}='worst', this would
%                                  correspond to the worst-case scenario
%                                  (Fig. 3a of the ICML paper)
%                                  If varargin{i+1}='best', this would
%                                  correspond to the best-case scenario
%                                  (Fig. 3b of the ICML paper)
%                                  (default: 'avg')
% 
% 
% OUTPUT
%  A figure diff.fig comparing the marginals and a figure speedup.fig
%  comparing the running times across the three methods and .mat file with 
% all the important output saved as: time_sp_*
% 
% 
% EXAMPLES
%
% Assume you want to randomly generate trees of varying size, where
% the cardinality among all nodes is 2.
% adabp_toy;
%
% Assume you want to randomly generate trees of varying size, where
% the cardinality among all nodes is 10.
% card = 10;
% adabp_toy(card);
% 
% Assume you want to randomly generate trees of varying size, where
% the cardinality among all nodes is 2 and AdaBP theoretically outperforms
% RCTreeBP by a ratio of 3.
% card = 2;
% adabp_toy(card,'speedup_rt',3);
% 
% Assume you want to randomly generate trees of varying size, where
% the cardinality among all nodes is 2 and the unbalanced tree ratio is 0.5
% card = 2;
% adabp_toy(card,'speedup_rt',NaN,'unbal_rt',0.5);
% 
% Assume you want to generate the worst-case scenario, where the 
% cardinality among all nodes is 2.
% card = 2;
% adabp_toy(2,'scenario','worst');
% Assume you want to generate the best-case scenario, where the 
% cardinality among all nodes is 2.
% card = 2;
% adabp_toy(2,'scenario','best');
% 
% Author: geopapa
% $ Date: 2015/06/29 18:21:33 $

    close all;

    N_arr     = 10.^[1  2  3  4];    % number of latent (hidden) nodes
    % N_arr     = 10.^[1  2  3];     % Note: it might take some time to run it on trees of 10^4 nodes.
    pmax      = [ 2   2   2     2];  % maximum number of children per node
    R_arr     = [10  10   10    2];  % number of different runs
    Nobs      = 50;                  % We assume that a maximum of Nobs observations (measurements) are drawn from each latent (hidden) node 
    disp_intv = 1;                   % The rate that we wish to have printed progress of the algorithm

    if nargin == 0,  card = 2;  end

    if any(strcmp(varargin,'speedup_rt'))
        speedup_rt = varargin{find(strcmp(varargin,'speedup_rt'))+1};
        if speedup_rt <= 0
            error('speedup_rt should be a positive number.');
        end
    else
        speedup_rt = NaN;
    end
    
    if any(strcmp(varargin,'unbal_rt'))
        unbal_rt = varargin{find(strcmp(varargin,'unbal_rt'))+1};
        if unbal_rt < 0 || unbal_rt > 1
            error('unbal_rt should be a number between 0 and 1.');
        end
    else
        unbal_rt = 0.7;
    end
    
    if any(strcmp(varargin,'scenario'))
        scenario = varargin{find(strcmp(varargin,'scenario'))+1};
        if ~(strcmp(scenario,'best') || strcmp(scenario,'worst'))
            error('scenario can only take value best or worst.');
        end
    else
        scenario = 'avg';
    end
    
    diff01 = zeros(max(R_arr),length(N_arr));
    diff12 = zeros(max(R_arr),length(N_arr));
    diff02 = zeros(max(R_arr),length(N_arr));
    time0  = nan(max(R_arr),length(N_arr));
    time1  = nan(max(R_arr),length(N_arr));
    time2  = nan(max(R_arr),length(N_arr));

    avg_path_len = nan(max(R_arr),length(N_arr));

    dateStart = datestr(now);
    t_s = tic;
    n = 1;
    for N = N_arr
        disp(['N(', num2str(n), ')=', num2str(N)]);
        for r = 1:R_arr(n)
            
            if strcmp(scenario,'avg')
                %% Construct the tree
                E     = zeros(N-1,2);
                Ne    = size(E,1);
                node  = 1;
                cnt   = 1;
                stack = [];
                while true
                    if rand < unbal_rt
                        num_suggested_children = 1;
                    else
                        num_suggested_children = pmax(n);
                    end            
                    numchild = min(num_suggested_children, Ne-cnt+1);

                    if numchild > 0
                        children = cnt:cnt+numchild-1;
                        children = children + 1;
                        children(children>N) = [];
                        E(cnt:cnt+numchild-1,:) = [repmat(node,numchild,1) children(:)];
                        cnt      = cnt + numchild;
                        stack    = [stack, children];
                        node     = stack(1);
                        stack(1) = [];
                    end

                    if cnt > Ne
                        break;
                    end
                end

                %% Assume we obtain one measurement per latent (hidden) node
                M = N;

                %% Determine the path lengths between all pairs of nodes
                pair_dist = eval_pair_dist(N, E);

                %% Construct the measurement order w
                if ~isnan(speedup_rt)
                    w = constructMeasOrder(N, M, card, speedup_rt, pair_dist);
                else
                    w = randperm(M);
                end
            elseif strcmp(scenario,'best')
                %% Construct the star graph
                E = [ones(1,N-1); 2:N];
                E = E';
                
                %% Assume we obtain one measurement per latent (hidden) node
                M = N;
                
                %% Construct the measurement order w
                w = randperm(M);
                
                %% Determine the path lengths between all pairs of nodes
                pair_dist = eval_pair_dist(N, E);
            else  % worst-case scenario
                %% Construct the Markov chain
                E = [(1:N-1); 2:N];
                E = E';
                
                %% Assume we obtain one measurement per latent (hidden) node
                M = N;
                
                %% Construct the measurement order w
                w = constructWorstCaseMeasOrder(N,M);
                
                %% Determine the path lengths between all pairs of nodes
                pair_dist = eval_pair_dist(N, E);
            end
            
            consec_dist = zeros(1,length(w)-1);
            for i = 1:length(consec_dist),  consec_dist(i) = pair_dist(w(i),w(i+1));  end
            disp(['    max(pair_dist) = ', num2str(max(pair_dist(:))), '  /  max(w1,w2) = ', num2str(max(consec_dist)), '  /  total # of msg = ', num2str(sum(consec_dist)), '  /  speedup_rt = ', num2str(speedup_rt), '  /  card = ', num2str(card)]);
            
            %% Set v to w (for convenience, this can be changed to an arbitrary sequence)
            v = w;
            
            %% Create the measurements
            y = randi(Nobs,1,M);

            %% Construct the potentials
            Ne      = size(E,1);
            phi     = cell(1,N);
            psi     = cell(1,Ne);
            phi_obs = cell(1,N);
            for i = 1:N,   phi{i}     = .005 + rand(card,1);     end
            for e = 1:Ne,  psi{e}     = .005 + rand(card,card);  end
            for i = 1:N,   phi_obs{i} = .005 + rand(card,Nobs);  end

            %% Choose a root randomly
            root = randi(N);

            %% Compare timing for standard BP
            disp('    BP running...');
            b0_v = cell(1,M);
            phi0 = phi;
            tic;
            for j = 1:M
                phi0{w(j)} = phi0{w(j)}.*phi_obs{w(j)}(:,y(j));
                b          = bp(E, phi0, psi, [], 'root', root);
                b0_v{j}    = b{v(j)};
            end
            time0(r,n) = toc;

            %% Compare timing for Adaptive Exact Inference by Sümer et al. (2011); RCTreeBP
            disp('    RCTreeBP running...');
            b1_v = cell(1,M);
            phi1 = phi;
            A    = sparse([1:N, E(:)'], [1:N, repmat(N+1:N+Ne,1,2)], 1);  % Adjacent matrix
            F    = [phi, psi];
            tic;
            % Build RC-tree structure for Sum-Product
            T = rctreeSP(A,F);
            for j = 1:M
                phi1{w(j)} = phi1{w(j)}.*phi_obs{w(j)}(:,y(j));
                T          = updateFactorSP(T, w(j), phi1{w(j)});
                b1_v{j}    = marginal(T,v(j));
            end
            time1(r,n) = toc;

            %% Compare timing for Adaptive BP; AdaBP
            disp('    AdaBP running...');
            b2_v = cell(size(b1_v));
            path_len_ww = zeros(1,M-1);
            path_len_wv = zeros(1,M);
            tic;
            abp = AdaBP(E,phi,psi,'root',root);
            for j = 1:M
                if j > 1
                    plw = abp.propagate(w(j-1), w(j), true);  % Propagate messages from w_{j-1} to w_j
                    path_len_ww(j-1) = plw;
                end
                abp.update(w(j),phi_obs{w(j)}(:,y(j)));  % Update the node potential at w_j

                plv = abp.propagate(w(j), v(j));  % Propagate messages from w_j to v_j
                path_len_wv(j) = plv;

                b2_v{j} = abp.eval_mrg(v(j));  % Evaluate the marginal at v_j
            end
            time2(r,n) = toc;

            % If v is same as w (we set that on purpose), we should not take it
            % into account for avegage path length estimation
            if all(path_len_wv==0),  path_len_wv = [];  end
            avg_path_len(r,n) = mean([path_len_ww,path_len_wv]);

            %% Compare the results of standard BP, RCTreeBP and AdaBP
            diff = 0;
            for j = 1:M
                diff = diff + norm(b0_v{j}-b1_v{j}, 1);
            end
            diff01(r,n) = mean(diff);

            diff = 0;
            for j = 1:M
                diff = diff + norm(b0_v{j}-b2_v{j}, 1);
            end
            diff02(r,n) = mean(diff);

            diff = 0;
            for j = 1:M
                diff = diff + norm(b1_v{j}-b2_v{j}, 1);
            end
            diff12(r,n) = mean(diff);

            if mod(r,disp_intv) == 0
                t_e = toc(t_s);
                disp(['    Sample # ', num2str(r), '/', num2str(R_arr(n)), ' has been processed. Start: ', dateStart, '. Elapsed time (', datestr(now), '): ', num2str(t_e), ' sec.']);
            end
            save(['time_sp_card_', num2str(card), '_speedup_', num2str(speedup_rt), '.mat'],'speedup_rt','card','N_arr','avg_path_len','diff01','diff02','diff12','time0','time1','time2');
        end
        n = n + 1;
    end

    if size(time0,1) == 1
        time0 = [time0; time0];
        time1 = [time1; time1];
        time2 = [time2; time2];
    end

    avg_path_len_mu       = nanmean(avg_path_len);
    rt_time1_to_time2_emp = nanmean(time1)./nanmean(time2);

    save(['time_sp_card_', num2str(card), '_speedup_', num2str(speedup_rt), '.mat'],'avg_path_len_mu','rt_time1_to_time2_emp','-append');

    disp('');
    disp(['speedup_rt = ', num2str(speedup_rt)]);
    disp(['card = ',  num2str(card)]);
    disp(['N = ',     num2str(N_arr)]);
    disp(['card*log10(N) = ', num2str(card), '*[', num2str(log10(N_arr)), '] = ', num2str(card*log10(N_arr))]);
    disp('');
    disp(['Empirical avg path len = [', num2str(avg_path_len_mu), ']']);
    disp('');
    disp(['Empirical speedup ratio time(RCTreeBP)/time(AdaBP) = [', num2str(rt_time1_to_time2_emp), ']']);
    disp('');

    h1 = figure(1);
    plot(N_arr,nanmean(diff01), 'r');
    hold on;
    plot(N_arr,nanmean(diff02), 'g');
    plot(N_arr, nanmean(diff12), 'b');
    axis tight;
    set(gca,'FontSize',30);
    xlabel('N');
    ylabel('$\mathbf{E}\left[\sum\limits_{i=1}^{M} \sum\limits_{x_i} |b_i^{(1)}(x_i) - b_i^{(2)}(x_i)|\right]$','Interpreter','LaTex');
    legend({'diff(BP,RCTreeBP)','diff(RCTreeBP,AdaBP)'}, 'Location', 'NorthWest','FontSize',26)
    saveas(h1,'diff.fig');

    h2 = figure(2);
    ratio0_2 = time0./time2;
    ratio1_2 = time1./time2;
    ratio0_2 = log10(ratio0_2);
    ratio1_2 = log10(ratio1_2);
    rt_mean  = [nanmean(ratio0_2); nanmean(ratio1_2)]';
    rt_std   = [nanstd(ratio0_2);  nanstd(ratio1_2)]';
    b        = bar(rt_mean);
    if ~verLessThan('matlab', '8.4')
        b(1).FaceColor = [211 211 211]/255;
        b(2).FaceColor = 'k';
    end
    hold on;
    [numDataPts, numComparisons] = size(rt_mean);
    if numDataPts > 1
        hErrbar = zeros(1,numComparisons);
        for j = 1:numComparisons
            % Extract the x location data needed for the errorbar plots:
            if verLessThan('matlab', '8.4'),  xLoc = get(get(b(j),'children'),'xdata');
            else                              xLoc = b(j).XData + [b(j).XOffset];
            end
            hErrbar(j) = errorbar(mean(xLoc,1), rt_mean(:,j), rt_std(:,j), rt_std(:,j), '.r');
            set(hErrbar(j), 'marker', 'none')
        end
    else
        % Extract the x location data needed for the errorbar plots:
        if verLessThan('matlab', '8.4'),  xLoc = get(get(handles.bar,'children'),'xdata');
        else                              xLoc = handles.bar.XData + [handles.bar.XOffset];
        end
        hErrbar = errorbar(mean(xLoc,1), rt_mean, rt_std, rt_std, '.r');
        set(hErrbar, 'marker', 'none')
    end
    axis tight;
    hold off;
    set(gca,'FontSize',30);
    xtickname = cellstr(num2str(round(log10(N_arr(:))), '10^%d'));
    set(gca,'XTickLabel',xtickname);
    max_y_tick = ceil(nanmax([nanmean(ratio0_2), nanmean(ratio1_2)]));
    min_y_tick = floor(nanmin([nanmean(ratio0_2), nanmean(ratio1_2)]));
    set(gca,'YLim', [min_y_tick max_y_tick]);
    set(gca,'YTick', min_y_tick:max_y_tick);
    yticks = round(get(gca,'YTick'));
    ytickname = cellstr(num2str(yticks(:), '10^{%d}'));
    set(gca,'YTickLabel',ytickname);

    if any(ratio1_2(:,end) < 0)  % RCTreeBP is faster than AdaBP
        text(1, -1, ['$|\mathcal{X}|=', num2str(card), '$'], 'interpreter','latex','FontSize',30);
    else
        text(1, 2, ['$|\mathcal{X}|=', num2str(card), '$'], 'interpreter','latex','FontSize',30);
    end
    
    xlabel('$N$', 'interpreter','latex');
    ylabel('Speedup ratio');
    legend({'t(BP)/t(AdaBP)', 't(RCTreeBP)/t(AdaBP)'},'FontSize',26, 'Location','NorthWest');
    %tightfig;
    saveas(h2,'speedup.fig');
end

% This function evaluates the pairwise distances between all (latent) nodes in the tree
function pair_dist = eval_pair_dist(N, E)
    pair_dist           = zeros(N);
    nghb                = GraphUtils.find_nghb(E);
    ch                  = GraphUtils.get_ch(nghb,1);
    [Eul, De, H, D] = GraphUtils.eulertour(1, ch);
    %treeplot(pa);
    Ii = rmq(De);
    for ii = 1:N-1
        for jj = ii+1:N
            if H(ii) < H(jj)
                low = H(ii);
                upp = H(jj);
            else
                low = H(jj);
                upp = H(ii);
            end

            kk = floor(log2(upp - low + 1));
            if De(Ii(low,kk+1)) <= De(Ii(upp-2^kk+1,kk+1))
                lca_idx = Ii(low,kk+1);
            else
                lca_idx = Ii(upp-2^kk+1,kk+1);
            end
            lca = Eul(lca_idx);

            pair_dist(ii,jj) = D(ii)+D(jj)-2*D(lca);
        end
    end
    pair_dist = pair_dist + pair_dist';
end

% This function constructs the measurement order w such that the speedup 
% ratio of AdaBP over RCTreeBP is approximately speedup_rt
function w = constructMeasOrder(N, M, card, speedup_rt, pair_dist)
    max_dist = round(card*log10(N)/speedup_rt);  
    w        = zeros(1,M);  % measurement order
    % Arbitrarily choose the first element
    w(1) = randperm(N,1);
    for j = 1:M-1
        dist_nxt_node = max_dist;
        Idist     = pair_dist(w(j),:) == dist_nxt_node;
        Idist_pl1 = pair_dist(w(j),:) == dist_nxt_node+1;
        Idist_mn1 = pair_dist(w(j),:) == dist_nxt_node-1;
        Idist     = Idist | Idist_pl1 | Idist_mn1;
        idx       = find(Idist);
        idx(randperm(length(idx))) = idx;
        if isempty(idx)
            [~, w(j+1)] = max(pair_dist(w(j),:));
        else
            w(j+1) = idx(1);
        end
    end
end

% This function constructs a worst-case scenario measurement order w, where
% the distace of consecutive elements in w is on the order of N.
function w = constructWorstCaseMeasOrder(N, M)
    w    = zeros(1,M);  % measurement order
    w(1) = randi(N);
    for j = 2:M
        if w(j-1) < N/2
            start_pos = min(w(j-1) + ceil(2*N/3),N);
            w(j)      = min(start_pos + randi(N-start_pos+1),N);
        else
            start_pos = max(w(j-1) - ceil(2*N/3),1);
            w(j)      = max(start_pos - randi(start_pos),1);
        end
    end
end