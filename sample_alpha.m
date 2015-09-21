function new_alpha = sample_alpha(alpha,theta,Component_Num,num_reg_used)
    
    alpha0 = Dirichlet_mle(theta, Component_Num);
    sum_log_theta = sum(log(theta));
    relax = 9;
    a = - sum_log_theta/relax;
    b = - alpha0'*relax./sum_log_theta;
    
    term = (sum_log_theta-1)*(alpha-1) + num_reg_used*( log(gamma(sum(alpha))) - sum(log(gamma(alpha))) );
    
    for kk = 1:Component_Num
        z = gamrnd(a(kk), b(kk)); % Gamma (independent) proposal
        alpha_kk = alpha(kk);
        alpha(kk) = z;
        term_z = (sum_log_theta-1)*(alpha-1) + num_reg_used*(log(gamma(sum(alpha))) - sum(log(gamma(alpha))) );
        A = exp(term_z - term + (a(kk)-1)*(log(alpha_kk)-log(z)) + (z - alpha_kk)/b(kk));
        u = rand(1);
        if u < A
            term = term_z;
        else
            alpha(kk) = alpha_kk;
        end
    end
    
    new_alpha = alpha;
    
end