function new_kappa = update_kappa(tau,i,j,num_reg_used,Method)
    
    % Update kappa
    tau_2d = sum((tau(i)-tau(j)).^2)/2; % each pair was computed twice in raw summation
    
    if tau_2d == 0
        tau_2d = 3e3;
    end
   
    if strcmp(Method,'MCMC')
        new_kappa = gamrnd((num_reg_used-1)/2, 2/tau_2d); % sample from kappa posterior
    else 
        new_kappa = (num_reg_used-3)/tau_2d; % kappa posterior mode
    end

end