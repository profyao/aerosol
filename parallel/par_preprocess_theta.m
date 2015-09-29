function [old_residp,old_thetap,taup,old_theta_neighbor,regp,smartp] = par_preprocess_theta(x,y,p,old_resid,old_theta,tau,i,j,reg,smart,kf,add_limit,const)

    xp = x(p);
    yp = y(p); 
    old_residp = old_resid(:,p);
    old_thetap = old_theta(:,p);
    taup = tau(p);

    id = find(j == p);
    if ~isempty(id)
        neighbor = i(id);
        old_theta_neighbor = old_theta(:,neighbor);
    else
        old_theta_neighbor = [];
    end

    [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf,add_limit);

end