function [atm_path,surface,resid] = get_resid(tau,theta,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit)
    
    [atm_path,~] = get_model(tau,theta,ExtCroSect,CompSSA,smartp,const,add_limit);

    %upbd = atm_path + surf_lim > regp.min_equ_ref;
    upbd = atm_path > 0.9 * regp.min_equ_ref;
    
    % simulation code
    % upbd = zeros(9,4);
    % simulation end
    
    if ~any(upbd(:)) && r > 1100
        [resid,surface] = get_resid_eof(regp,atm_path,const);
        % simulation code
        %resid = regp.mean_equ_ref - atm_path;
        %resid = resid(:);
        %surface = zeros(const.NChannel,1);
        % simulatino end
    elseif ~any(upbd(:)) && r == 1100
        [resid,surface] = get_resid_eof2(regp,atm_path,const);
    else
        resid = Inf * ones(const.NChannel,1);
        surface = NaN * ones(const.NChannel,1);
    end
    
    atm_path = atm_path(:);

end