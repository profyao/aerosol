function [atm_path,surf,resid] = get_resid(tau,theta,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit)
    
    [atm_path,~] = get_model(tau,theta,ExtCroSect,CompSSA,smartp,const,add_limit);

    %upbd = atm_path + surf_lim > regp.min_equ_ref;
    upbd = atm_path > regp.min_equ_ref * 0.9;
    
    if ~any(upbd(:)) && r > 1100
        [resid,surf] = get_resid_eof(regp,atm_path,const);
    elseif ~any(upbd(:)) && r == 1100
        [resid,surf] = get_resid_eof2(regp,atm_path,const);
    else
        resid = Inf * ones(const.NChannel,1);
        surf = NaN * ones(const.NChannel,1);
    end
    
    atm_path = atm_path(:);

end