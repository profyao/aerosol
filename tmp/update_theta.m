function [new_theta, new_resid] = update_theta(current, i, j, x, y, smart, reg, ExtCroSect, CompSSA)

    new_theta = current.theta;
    new_resid = current.resid;
    
    for p = 1:reg.num_reg_used
        
        neighbor = find(j == p);
        
        if ~isempty(neighbor)
            mu = mean(current.theta(:,i(neighbor)),2);
        else
            mu = mean(current.theta,2);
        end

        theta = gamrnd(mu,1);
        theta = theta / sum(theta);
        
        [~,~,resid] = get_resid(current.tau(p),theta,x(p),y(p),reg,smart,ExtCroSect,CompSSA);

        if nansum(current.resid(:,p).^2./current.sigmasq) > nansum(resid.^2./current.sigmasq) || ...
                ( isinf(current.resid(1,p)) && isinf(resid(1)) )
            new_theta(:,p) = theta;
            new_resid(:,p) = resid;
        end
    end

end
