function [new_theta, new_resid] = par_update_theta(old_theta,tau,old_resid,sigmasq,alpha,...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, kf, par, add_limit,const)
    
    if par == true
        
        new_theta = old_theta;
        new_resid = old_resid;
    
        parfor p = 1:reg.num_reg_used

            [old_residp,old_thetap,taup,old_theta_neighbor,regp,smartp] = par_preprocess_theta(x,y,p,old_resid,old_theta,tau,i,j,reg,smart,kf,const);
            
            [new_theta(:,p),new_resid(:,p)] = par_update_theta_pixel(old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
                regp,smartp,ExtCroSect,CompSSA,kf,add_limit,const);   
        end
        
    else
    
        for p = 1:reg.num_reg_used

            [old_residp,old_thetap,taup,old_theta_neighbor,regp,smartp] = par_preprocess_theta(x,y,p,old_resid,old_theta,tau,i,j,reg,smart,kf,const);

            [old_theta(:,p),old_resid(:,p)] = par_update_theta_pixel(old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
                regp,smartp,ExtCroSect,CompSSA,kf,add_limit,const);   
        end
        
        new_theta = old_theta;
        new_resid = old_resid;
             
    end

end
