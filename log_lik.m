function [y,num] = log_lik(current,i,j,Channel_Used,Method)
    
    succ = ~isinf(current.resid(1,:));
    P = length(unique(j));
    num = sum(succ);
    id = current.sigmasq~=0 & Channel_Used;
    
    if strcmp(Method,'CD-random-noprior') || strcmp(Method,'CD-noprior')
        y = - sum(0.5*log(current.sigmasq(id)).* sum(~isnan(current.resid(id,succ)),2)+2) - ...
            0.5 * nansum(nansum(current.resid(id,succ).^2./repmat(current.sigmasq(id),1,num)));
    else
        sum_log_theta = sum(log(current.theta'));
       
        y = - sum(0.5*log(current.sigmasq(id)).* sum(~isnan(current.resid(id,succ)),2)+2) - ...
        0.5 * nansum(nansum(current.resid(id,succ).^2./repmat(current.sigmasq(id),1,num))) - ...
        0.25 * current.kappa * sum((current.tau(i) - current.tau(j)).^2) ... 
        + (P-3)/2 * log(current.kappa) + ...
        (sum_log_theta-1)*(current.alpha-1) + P *( log(gamma(sum(current.alpha))) - sum(log(gamma(current.alpha))));
       
    end
    
end