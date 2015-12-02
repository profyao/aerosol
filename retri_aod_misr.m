function [chisq,chisq_misr,tau_misr,theta_misr] = retri_aod_misr(theta_misr_grid,r,sample,x,y,const,ExtCroSect,CompSSA,smart,reg)

    chisq = zeros(reg.num_reg_used,1);
    %term_kappa = zeros(reg.num_reg_used,1);
    chisq_misr = zeros(reg.num_reg_used,1);
    %term_kappa_misr = zeros(reg.num_reg_used,1);
    tau_misr = zeros(reg.num_reg_used,1);
    theta_misr = zeros(const.Component_Num,reg.num_reg_used);
    sigmasq = sample.sigmasq(:,end);
    pid = [119,255,240,237,267];
    theta = sample.theta;
    tau = sample.tau;

    parfor p=1:reg.num_reg_used
        thetap = theta(:,p);
        taup = tau(p);
        xp = x(p);
        yp = y(p);
        chisq(p) = 0.5 * get_chisq(r,taup,thetap,sigmasq,xp,yp,const,ExtCroSect,CompSSA,smart,reg);
        %term_kappa(p) = 0.5 * get_term_kappa(i1d,j1d,p,sample.tau,taup);
        [chisq_misr(p),tau_misr(p),theta_misr(:,p)] = get_misr_retri(r,theta_misr_grid,sigmasq,xp,yp,ExtCroSect,CompSSA,smart,reg,const);
        %term_kappa_misr(p) = 0.5 * get_term_kappa(i1d,j1d,p,sample.tau,tau_misr(p));
    end

end