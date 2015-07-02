function [mu1, var1, mu2, var2] = adabp_loopy_toy()
    N = 18;
    M = 1000;
    
    [h, J] = form_h_J();

    A = [4 9 11 13 14 15];  % anchor nodes
    A = A(:)';
    F = [16 17 18];         % FVS
    T = setdiff(1:N,F);     % Tree part of graph
    
    
    nghb_F = {[4  9 ]; 
              [11 14]; 
              [13 15]};

    nghb = nan(N,1);  nghb(4) = F(1);  nghb(9) = F(1);  nghb(11) = F(2);  nghb(14) = F(2);  nghb(13) = F(3);  nghb(15) = F(3);

    isT    = true(1,N);
    isT(F) = false;

    h_T  = h(T);
    J_TT = J(T,T);
    h_F  = h(F);
    J_FF = J(F,F);

    w    = randi(N,1,M);
    v    = randi(N,1,M);
    Y    = randn(1,M);
    C    = randi(3,1,M);
    muW  = zeros(1,M);
    R    = randi(4,1,M);
    mu1  = nan(1,M);
    var1 = nan(1,M);
    h_cp = h;
    J_cp = J;
    for j = 1:M
        % Incorporate the measurement
        [h_cp(w(j)), J_cp(w(j),w(j))] = set_h_J(h_cp(w(j)), J_cp(w(j),w(j)), j, Y, C, muW, R);
        S       = J_cp\eye(size(J_cp));
        mu      = S*h_cp;
        mu1(j)  = mu(v(j));
        var1(j) = S(v(j),v(j));
    end
    
    exclIdx    = randperm(M,600);
    v(exclIdx) = 0;
    
    mu2      = nan(1,M);
    var2     = nan(1,M);
    abp_T    = AdaBP(h_T, J_TT);
    abp_p{1} = AdaBP(J(T,F(1)), J_TT);
    abp_p{2} = AdaBP(J(T,F(2)), J_TT);
    abp_p{3} = AdaBP(J(T,F(3)), J_TT);
    
    % Evaluate the partial means and ``feedback gains''
    [mu_A, g_A] = eval_mu_A_g_A(abp_T, abp_p, A, F, N);
    
    % Determine mu_F, S_F
    [mu_F, S_F] = eval_mu_F_S_F(h_F, J_FF, J, F, nghb_F, mu_A, g_A);
    
    % Start the incremental process
    prev_w = 0;
    for j = 1:M
        if j > 1
            %% Ensure consistency by propagating to the next measurement node
            if isT(w(j-1))  % Send messages from w_{j-1} to w_{j} in T
                prev_w = w(j-1);
            end
            % else: Send messages from previous w belonging in T to w_{j} in T
            
            if isT(w(j)) && prev_w~=0
                for f = 1:length(F)
                    abp_p{f}.propagate(prev_w, w(j), true);  % h_{i->j}^p, J_{i->j}^T
                end
                abp_T.propagate(prev_w, w(j), true);  % h_{i->j}^T, J_{i->j}^T
            end
        end
        
        %% Incorporate the measurement at w_{j}
        [h(w(j)), J(w(j),w(j))] = set_h_J(h(w(j)), J(w(j),w(j)), j, Y, C, muW, R);
        if isT(w(j))  % If measurement is from T, update h_{w(j)}, J_{w(j),w(j)}
            for f = 1:length(F)
                abp_p{f}.setNodePot(w(j), [], J(w(j),w(j)));
            end
            abp_T.setNodePot(w(j), h(w(j)), J(w(j),w(j)));
        else  % If measurement is from F (FVS), update the potential vector and information matrix
            h_F  = h(F);
            J_FF = J(F,F);
        end
        
        if v(j) ~= 0  % It means that there is a marginal to evaluate at this stage
            %% Propagate to marginal node v_{j}
            if isT(w(j))
                cur_w = w(j);    % Send messages from w_{j}      belonging in T to v_{j} in T
            else
                cur_w = prev_w;  % Send messages from previous w belonging in T to v_{j} in T
            end

            % Send messages from cur_w to all anchor nodes
            if cur_w ~= 0
                for f = 1:length(F)
                    abp_p{f}.propagate(cur_w, A);  % h_{i->j}^p, J_{i->j}^T
                end
                abp_T.propagate(cur_w, A);         % h_{i->j}^T, J_{i->j}^T
            end

            % Evaluate the partial means and ``feedback gains''
            [mu_A, g_A] = eval_mu_A_g_A(abp_T, abp_p, A, F, N);

            % Determine mu_F, S_F
            [mu_F, S_F] = eval_mu_F_S_F(h_F, J_FF, J, F, nghb_F, mu_A, g_A);

            % Introduce the effect of the ``fictitious'' revised potential vectors in A
            if isT(v(j))
                % Evaluate h_tilde and propagate the h messages
                h_A_tilde = eval_h_A_tilde(h, A, J, mu_F, nghb, F);
                h_A_old   = h(A);
                cnt = 1;
                for anchor = A
                    abp_T.setNodePot(anchor, h_A_tilde(cnt));
                    cnt = cnt + 1;
                end
                abp_T.propagate(A, v(j));
                
                for f = 1:length(F)
                    abp_p{f}.propagate(cur_w, v(j));  % h_{i->j}^p, J_{i->j}^T
                end
            end

            %% Evaluate the marginal of the marginal node
            if isT(v(j))
                [h_mrg,J_mrg] = abp_T.eval_mrg(v(j));  % Messages J_{i->j}^T are unaffected from the change in revised potential vectors h_tilde
                g_v           = zeros(length(F),1);
                for f = 1:length(F)
                    [h_tmp, J_tmp] = abp_p{f}.eval_mrg(v(j));
                    g_v(f)         = J_tmp\h_tmp;
                end
                var2(j) = J_mrg\eye(size(J_mrg)) + g_v'*S_F*g_v;
                mu2(j)  = J_mrg\h_mrg;
            else
                var2(j) = S_F(v(j)-F(1)+1,v(j)-F(1)+1);
                mu2(j)  = mu_F(v(j)-F(1)+1);
            end

            % Scratch the effect of the revised potential vectors in A, since these were fictitious and not real changes
            if isT(v(j))
                cnt = 1;
                for anchor = A
                    abp_T.setNodePot(anchor, h_A_old(cnt));
                    cnt = cnt + 1;
                end
                abp_T.propagate(A, v(j));
            end
        end
    end
        [min(w), max(w)]
        [min(v(v>0)), max(v)]
    
    figure(1);
    subplot(2,1,1);
    plot(setdiff(1:M,exclIdx), mu1(setdiff(1:M,exclIdx)), 'g', setdiff(1:M,exclIdx), mu2(setdiff(1:M,exclIdx)), 'b', 'LineWidth', 1.5);
    set(gca, 'FontSize', 30);
    title('\mu');
    legend({'J^{-1}h', 'IBP'}, 'FontSize', 26);    
    axis tight;
    subplot(2,1,2);
    plot(setdiff(1:M,exclIdx), var1(setdiff(1:M,exclIdx)), 'g', setdiff(1:M,exclIdx), var2(setdiff(1:M,exclIdx)), 'b', 'LineWidth', 1.5);
    set(gca, 'FontSize', 30);
    title('var');
    legend({'J^{-1}', 'adaBP'}, 'FontSize', 26);    
    axis tight;
    
    figure(2);
    subplot(2,1,1);
    mu12_diff_abs  = abs(mu1-mu2);
    var12_diff_abs = abs(var1-var2);
    stem(setdiff(1:M,exclIdx), mu12_diff_abs(setdiff(1:M,exclIdx)), '.r', 'Marker', 'None', 'LineWidth', 1.5);
    set(gca, 'FontSize', 30);
    title('$|\mu_1-\mu_2|$','Interpreter','LaTex');
    axis tight;    
    subplot(2,1,2);
    stem(setdiff(1:M,exclIdx), var12_diff_abs(setdiff(1:M,exclIdx)), '.r', 'Marker', 'None', 'LineWidth', 1.5);
    set(gca, 'FontSize', 30);
    title('$|\sigma_1^2-\sigma_2^2|$','Interpreter','LaTex');
    xlabel('# of iterations');
    axis tight;
end

function [h, J] = form_h_J()
    h = ones(N,1);
    J = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
         1  2  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;
         0  0  1  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0;
         0  0  1  5  0  0  0  1  0  0  0  0  0  0  0  1  0  0;
         0  0  0  0  4  1  0  0  0  0  0  0  0  0  0  0  0  0;
         0  1  0  0  1  1  1  0  0  0  1  0  0  0  0  0  0  0;
         0  0  1  0  0  1  3  0  0  0  0  0  0  0  0  0  0  0;
         0  0  0  1  0  0  0  6  1  0  0  0  0  0  0  0  0  0;
         0  0  0  0  0  0  0  1  3  0  0  0  0  0  0  1  0  0;
         0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0;
         0  0  0  0  0  1  0  0  0  1  1  1  0  0  0  0  1  0;
         0  0  0  0  0  0  0  0  0  0  1  2  1  1  0  0  0  0;
         0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  0  0  1;
         0  0  0  0  0  0  0  0  0  0  0  1  0  3  0  0  1  0;
         0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  1;
         0  0  0  1  0  0  0  0  1  0  0  0  0  0  0  1  0  0;
         0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  3  0;
         0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  2];
    J(1,2)   =  0.1;  J(2,1)   =  0.1;
    J(13,18) = -0.6;  J(18,13) = -0.6;
    J(3,4)   =  0.2;  J(4,3)   =  0.2;
    J(10,11) = -0.1;  J(11,10) = -0.1;
    J(13,15) =  0.7;  J(15,13) =  0.7;
    J(11,17) = -0.3;  J(17,11) = -0.3;
    J(4,8)   =  0.7;  J(8,4)   =  0.7;
    J(5,6)   = -0.2;  J(6,5)   = -0.2;
    J(6,7)   =  0.1;  J(7,6)   =  0.1;
    J(3,7)   = -0.5;  J(7,3)   = -0.5;
    J(2,6)   =  0.1;  J(6,2)   =  0.1;
    J(8,9)   =  0.2;  J(9,8)   =  0.2;
    J(4,16)  = -0.4;  J(16,4)  = -0.4;
    J(15,18) = -0.1;  J(18,15) = -0.1;
    J(11,12) = -0.1;  J(12,11) = -0.1;
    J(12,13) =  0.1;  J(13,12) =  0.1;
    J(6,11)  =  0.1;  J(11,6)  =  0.1;
end

function [mu_A, g_A] = eval_mu_A_g_A(abp_T, abp_p, A, F, N)
    mu_A = nan(N,1);
    g_A  = nan(N,length(F));
    for anchor = A
        [h_tmp, J_tmp] = abp_T.eval_mrg(anchor);  % J_tmp is the same across abp_T and all abp_p
        mu_A(anchor)   = J_tmp\h_tmp;
        for f = 1:length(F)
            h_tmp          = abp_p{f}.eval_mrg(anchor);
            g_A(anchor, f) = J_tmp\h_tmp;
        end
    end
end

function [mu_F, S_F] = eval_mu_F_S_F(h_F, J_FF, J, F, nghb_F, mu_A, g_A)
    h_F_hat = h_F;
    J_F_hat = J_FF;
    for p = 1:length(F)
        for q = 1:length(F)
            J_F_hat(p,q) = J_F_hat(p,q) - sum( J(nghb_F{p},F(p)).*g_A(nghb_F{p},q) );
        end
        h_F_hat(p) = h_F_hat(p) - sum( J(nghb_F{p},F(p)).*mu_A(nghb_F{p}) );
    end
	S_F  = J_F_hat\eye(size(J_F_hat));
	mu_F = S_F*h_F_hat;
end

function h_A_tilde = eval_h_A_tilde(h, A, J, mu_F, nghb, F)
    h_A_tilde = zeros(length(A),1);
    cnt = 1;
    for anchor = A
        h_A_tilde(cnt) = h(anchor) - J(anchor,nghb(anchor))*mu_F(nghb(anchor)-F(1)+1);
        cnt = cnt + 1;
    end
end

function [h, J] = set_h_J(h, J, j, Y, C, muW, R)
    tmp1 = C(j)'/R(j);
    tmp2 = tmp1*C(j);
    tmp2 = (tmp2 + tmp2')/2;
    h    = h + tmp1*(Y(j)-muW(j));
    J    = J + tmp2;
    J    = (J+J')/2;
end
