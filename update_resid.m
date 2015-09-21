function [atm_path,surf,resid] = update_resid(current, x, y, smart, reg, ExtCroSect, CompSSA, kf, add_limit, const)
    
    atm_path = NaN*ones(const.NChannel,reg.num_reg_used);
    surf = NaN*ones(const.NChannel,reg.num_reg_used);
    resid = NaN*ones(const.NChannel,reg.num_reg_used);
    
    for p = 1:reg.num_reg_used
        
        xp = x(p);
        yp = y(p); 
        thetap = current.theta(:,p); 
        taup = current.tau(p);
        
        [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf);
        
        [atm_path(:,p),surf(:,p),resid(:,p)] = get_resid(taup,thetap,xp,yp,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
    end

end