function sample = par_aod_retri_search(x,y,reg,smart,CompModNum,MixSSA,CompFrac,ExtCroSect,const,add_limit,kf)

    sample.tau = NaN*ones(reg.num_reg_used,1);
    sample.theta = NaN*ones(3,reg.num_reg_used);
    
    for p = 1:reg.reg_num_used
        
        xp = x(p);
        yp = y(p);
        resid = NaN*ones(const.Model_OpticalDepthFinerGridLen,const.Model_MixtureDim);
        
        for mixture = 1:const.Model_MixtureDim
            
            const.Component_Particle = int8(CompModNum(:, mixture));
            const.Component_Num = length(const.Component_Particle);
            [regp,smartp] = extract_pixel(xp,yp,reg,smart,const);
            
            for tau_id = 1:const.Model_OpticalDepthFinerGridLen
                tau = const.Model_OpticalDepthFinerGrid(tau_id);
                [~,~,resid(tau_id,mixture)] = get_resid_mixture(tau,mixture,regp,smartp,MixSSA,CompFrac,ExtCroSect,const,add_limit,kf);
                
            end
                
        end
        
    end

end