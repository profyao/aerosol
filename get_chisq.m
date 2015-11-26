function chisq = get_chisq(r,taup,theta,p,sigmasq,x,y,const,ExtCroSect,CompSSA,smart,reg,reg_sim)
    
    xp = x(p);
    yp = y(p);

    thetap = theta(:,p);
    
    reg.mean_equ_ref = reg_sim;
    reg.min_equ_ref = reg_sim;
    [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,r,0);
    [~,~,resid] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,r,0);
    
    chisq = nansum(resid.^2./sigmasq);
end