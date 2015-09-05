function new_sigmasq = update_sigmasq(resid,Method) 

    succ = ~isinf(resid(1,:));
    
    if strcmp(Method,'MCMC')
        new_sigmasq = nansum(resid(:,succ).^2,2)./chi2rnd(sum(~isnan(resid(:,succ)),2));
    else
        new_sigmasq = nansum(resid(:,succ).^2,2)./(sum(~isnan(resid(:,succ)),2)+2);
    end

end