function [atm_path,surf,resid] = update_resid(current, x, y, smart, reg, ExtCroSect, CompSSA, const)
    
    NChannel = const.NChannel;
    Band_Dim = const.Band_Dim;
    Cam_Dim = const.Cam_Dim;
    
    atm_path = NaN*ones(NChannel,reg.num_reg_used);
    surf = NaN*ones(NChannel,reg.num_reg_used);
    resid = NaN*ones(NChannel,reg.num_reg_used);
    
    parfor p = 1:reg.num_reg_used
        
        xp = x(p);
        yp = y(p); 
        thetap = current.theta(:,p); 
        taup = current.tau(p);
        
        channel_is_used = reshape(reg.channel_is_used(xp,yp,:,:),Band_Dim,Cam_Dim)';
        min_equ_ref = reshape(reg.min_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)';
        mean_equ_ref = reshape(reg.mean_equ_ref(xp,yp,:,:),Band_Dim,Cam_Dim)';
        eof = reshape(reg.eof(xp,yp,:,:,:),Band_Dim,Cam_Dim,Cam_Dim);
        max_usable_eof = reshape(reg.max_usable_eof(xp,yp,:),Band_Dim,1);
        
        [atm_path(:,p),surf(:,p),resid(:,p)] = get_resid(taup,thetap,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
            smart,ExtCroSect,CompSSA,const);
    end

end