function [atm_path,surf,resid] = get_resid_mixture(tau,mixture,regp,smartp,MixSSA,CompFrac,ExtCroSect,const,add_limit,kf)

    [atm_path,surf_lim] = get_model_mixture(tau,mixture,MixSSA,CompFrac,ExtCroSect,smartp,const,add_limit);

    upbd = atm_path + surf_lim > regp.min_equ_ref;
   
    if ~any(upbd(:)) && kf == false
        [resid,surf] = get_resid_eof(regp,atm_path,const);
    elseif ~any(upbd(:)) && kf == true
        [resid,surf] = get_resid_eof2(regp,atm_path,const);
    else
        resid = Inf * ones(const.NChannel,1);
        surf = NaN * ones(const.NChannel,1);
    end
    
    atm_path = atm_path(:);

end