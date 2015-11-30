function [chisq_misrp,tau_misrp,theta_misrp] = get_misr_retri(r,theta_misr_grid,sigmasq,xp,yp,ExtCroSect,CompSSA,smart,reg,const)
    chisq_grid = gen_tt_grid(r,theta_misr_grid,sigmasq,xp,yp,ExtCroSect,CompSSA,smart,reg,const);
    [chisq_misrp,ind] = min(chisq_grid(:));
    [tau_id,theta_id] = ind2sub([13,21],ind);
    tau_misrp = const.Model_OpticalDepthGrid(tau_id);
    theta_misrp = theta_misr_grid(:,theta_id);

end