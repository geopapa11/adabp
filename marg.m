function b = marg(E, phi, psi, v)
    card = cellfun(@length,phi);  % cardinality (alphabet) for each node
    N    = length(phi);
    K    = prod(card);
    seq  = zeros(K,N);
    for i = 1:N
        if i < N
            seq_i = repmat(1:card(i), prod(card(i+1:end)), 1);
        else
            seq_i = 1:card(i);
        end
        seq_i    = seq_i(:);
        seq_i    = repmat(seq_i, 1, prod(card(1:i-1)));
        seq_i    = seq_i(:);
        seq(:,i) = seq_i;
    end
    
    prob = zeros(K,1);
    for k = 1:K
        prob(k) = exp(loglik(seq(k,:), E, phi, psi));
    end
    
    b = zeros(card(v),1);
    for j = 1:card(v)
        b(j) = sum(prob(seq(:,v)==j));
    end
    b = b/sum(b);
end