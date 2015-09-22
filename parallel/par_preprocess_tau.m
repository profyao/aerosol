function [old_residp,old_taup,thetap,old_tau_neighbor,regp,smartp]...
    = par_preprocess_tau(x,y,p,old_resid,theta,i,j,old_tau,reg,smart,kf,const)

        xp = x(p);
        yp = y(p);
        old_residp = old_resid(:,p);
        old_taup = old_tau(p);
        thetap = theta(:,p);
        
        id = find(j == p);
        if ~isempty(id)
            neighbor = i(id);
            old_tau_neighbor = old_tau(neighbor);
        else
            old_tau_neighbor = [];
        end
        
        [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf);
 

end