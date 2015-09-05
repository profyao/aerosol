function [atm_path,surf,resid] = get_resid(tau,theta,x,y,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
            smart,ExtCroSect,CompSSA,const)

    NChannel = const.NChannel;
    
    atm_path = get_model(tau,theta,x,y,ExtCroSect,CompSSA,smart,const);

    upbd = atm_path > min_equ_ref;
   
    if ~any(upbd(:))
        [resid,surf] = get_resid_eof(channel_is_used,mean_equ_ref,eof,max_usable_eof,atm_path,const);
    else
        resid = Inf * ones(NChannel,1);
        surf = NaN * ones(NChannel,1);
    end
    
    atm_path = atm_path(:);

end