function [new_tau,new_resid] = par_update_tau(old_tau,theta,old_resid,kappa,sigmasq,delta, ...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, const)
    
    Band_Dim = const.Band_Dim;
    Cam_Dim = const.Cam_Dim;
    
    new_tau = old_tau;
    new_resid = old_resid;
        
    parfor p = 1:reg.num_reg_used
                
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
        
        channel_is_used = reshape(reg.channel_is_used(xp,yp,:,:),Band_Dim,Cam_Dim)'; %9x4
        min_equ_ref = reshape(reg.min_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)'; %9x4
        mean_equ_ref = reshape(reg.mean_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)'; %9x4
        eof = reshape(reg.eof(xp,yp,:,:,:),Band_Dim,Cam_Dim,Cam_Dim); %4x9x9
        max_usable_eof = reshape(reg.max_usable_eof(xp,yp,:),Band_Dim,1); %4x1
        
        [new_tau(p),new_resid(:,p)] = par_update_tau_pixel(xp,yp,old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
    channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,ExtCroSect,CompSSA,smart,const);        
        
    end
    
end