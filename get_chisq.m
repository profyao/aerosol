function chisq = get_chisq(r,taup,thetap,sigmasq,xp,yp,const,ExtCroSect,CompSSA,smart,reg)
    
    %load simulation data
    %reg.mean_equ_ref = reg_sim;
    %reg.min_equ_ref = reg_sim;
    
    [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,r,0);
    [~,~,resid] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,r,0);
    
    chisq = nansum(resid.^2./sigmasq);
end