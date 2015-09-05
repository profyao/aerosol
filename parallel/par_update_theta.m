function [new_theta, new_resid] = par_update_theta(old_theta,tau,old_resid,sigmasq,alpha,...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, const)

    Band_Dim = const.Band_Dim;
    Cam_Dim = const.Cam_Dim;
    
    new_theta = old_theta;
    new_resid = old_resid;
    
    parfor p = 1:reg.num_reg_used
        
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
        
        channel_is_used = reshape(reg.channel_is_used(xp,yp,:,:),Band_Dim,Cam_Dim)';
        min_equ_ref = reshape(reg.min_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)';
        mean_equ_ref = reshape(reg.mean_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)';
        eof = reshape(reg.eof(xp,yp,:,:,:),Band_Dim,Cam_Dim,Cam_Dim);
        max_usable_eof = reshape(reg.max_usable_eof(xp,yp,:),Band_Dim,1);
        
        [new_theta(:,p),new_resid(:,p)] = par_update_theta_pixel(xp,yp,old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
            channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,smart,ExtCroSect,CompSSA,const);
        
    end

end
