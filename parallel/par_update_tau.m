function [new_tau,new_resid] = par_update_tau(old_tau,theta,old_resid,kappa,sigmasq,delta, ...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, kf, par, add_limit, const)
    
    if par == true
        
        new_tau = old_tau;
        new_resid = old_resid;
        
        parfor p = 1:reg.num_reg_used

            [old_residp,old_taup,thetap,old_tau_neighbor,regp,smartp] = par_preprocess_tau(x,y,p,old_resid,theta,i,j,old_tau,reg,smart,kf,const);

            [new_tau(p),new_resid(:,p)] = par_update_tau_pixel(old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
            regp,smartp,ExtCroSect,CompSSA,kf,add_limit,const);        

        end
        
    else
        
        for p = 1:reg.num_reg_used

            [old_residp,old_taup,thetap,old_tau_neighbor,regp,smartp] = par_preprocess_tau(x,y,p,old_resid,theta,i,j,old_tau,reg,smart,kf,const);

            [old_tau(p),old_resid(:,p)] = par_update_tau_pixel(old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
            regp,smartp,ExtCroSect,CompSSA,kf,add_limit,const);        

        end
        
        new_tau = old_tau;
        new_resid = old_resid;
    end
        
    
end