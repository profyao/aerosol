function [new_tau,new_resid] = update_tau(current, delta, i, j, x, y, smart, reg, ExtCroSect, CompSSA)
    
    new_tau = current.tau;
    new_resid = current.resid;
    
    for p = 1:reg.num_reg_used
        
        if isinf(current.resid(1,p))
            tau = current.tau(p) * 0.8;
            smooth0 = 0;
            smooth1 = 0;
        else       
            neighbor = find(j == p);
            
            if ~isempty(neighbor)
                mu = 0.5 * (mean(current.tau(i(neighbor))) + current.tau(p));
                tau = mu + delta * randn(1);
                tau(tau<=1e-3)=1e-3;
                tau(tau>=3)=3;
                smooth0 = current.kappa * sum(current.tau(p) - current.tau(i(neighbor))).^2;
                smooth1 = current.kappa * sum(tau-current.tau(i(neighbor))).^2;
            else
                mu = 0.5 * (mean(current.tau) + current.tau(p));
                tau = mu + delta * randn(1);
                tau(tau<=1e-3)=1e-3;
                tau(tau>=3)=3;
                %smooth0 = kappa * (current.tau(p) - mean(current.tau))^2;
                %smooth1 = kappa * (tau-mean(current.tau))^2;
                smooth0=0;smooth1=0;
            end       
        end
        
        [~,~,resid] = get_resid(tau,current.theta(:,p),x(p),y(p),reg,smart,ExtCroSect,CompSSA);
        
        if nansum(current.resid(:,p).^2./current.sigmasq) + smooth0 > nansum(resid.^2./current.sigmasq) + smooth1 ||...
                ( isinf(current.resid(1,p)) && isinf(resid(1)) )
            new_tau(p) = tau;
            new_resid(:,p) = resid;
        end

    end
    
end