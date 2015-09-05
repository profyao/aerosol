function [y,num] = log_lik(current,i,j)
    
    succ = ~isinf(current.resid(1,:));
    P = length(unique(j));
    num = sum(succ);
    id = current.sigmasq~=0;
    y = - sum(0.5*log(current.sigmasq(id)).* sum(~isnan(current.resid(id,succ)),2)+2) - ...
        0.5 * nansum(nansum(current.resid(id,succ).^2./repmat(current.sigmasq(id),1,num))) - ...
        0.25 * current.kappa * sum((current.tau(i) - current.tau(j)).^2) + ...
        (P-3)/2 * log(current.kappa);
    
end