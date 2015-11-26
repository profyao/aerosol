function new_alpha = update_alpha(old_alpha,theta,Component_Num,num_reg_used)
    
    alpha = Dirichlet_mle(theta, Component_Num);
    sum_log_theta = sum(log(theta));
   
    l1 = (sum_log_theta-1)*(alpha-1) + num_reg_used*( log(gamma(sum(alpha))) - sum(log(gamma(alpha))) );
    l0 = (sum_log_theta-1)*(old_alpha-1) + num_reg_used*( log(gamma(sum(old_alpha))) - sum(log(gamma(old_alpha))) );
    
    if l1 > l0    
        new_alpha = alpha;
    else
        new_alpha = old_alpha;
    end
    
end