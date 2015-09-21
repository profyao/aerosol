function [atm_path,surf_lim] = get_model_mixture(tau,mixture,MixSSA,CompFrac,ExtCroSect,smartp,const,add_limit)

    atm_path = zeros(const.Cam_Dim,const.Band_Dim);
    surf_lim = zeros(const.Cam.Dim,const.Band.Dim);
    
    for band = 1:const.Band_Dim
        
        scale_factor = ExtCroSect(band, const.Component_Particle) ./ ExtCroSect(const.Band_Green, const.Component_Particle) * CompFrac(const.Band.Green,:,mixture)';
        fraction_band = reshape(CompFrac(band,:,mixture),1,const.Component_Num);

        ssa_v = CompSSA(band, const.Component_Particle); %1x3
        ssa_mixture = MixSSA(band,mixture);
        
        tau_band = double(tau*scale_factor);
        ss_band = smartp.tau_cam_ss(:,:,:,band);
        ms_band = smartp.tau_cam_ms(:,:,:,band);
        rayleigh = reshape(ms_band(1,:,:,band),const.Cam_Dim,const.Component_Num);
        
        if tau_band==0
            ss = reshape(ss_band(1,:,:),const.Cam_Dim,const.Component_Num);
            ms = rayleigh;
        else          
            ss = interpol_3d_mex(tau_band, ss_band); %9x8
            ms = interpol_3d_mex(tau_band, ms_band);
        end

        atm_path(:, band) = (rayleigh + ss) * fraction_band' +  (ms - rayleigh) * (ssa_mixture ./ ssa_v .* exp(- tau_band * ...
            abs(ssa_mixture - ssa_v)) .* fraction_band)'; %9x1
                
        if add_limit == true
            tau_ediff_band = smartp.tau_ediff(:,:,band);
            tau_tdiff_band = smartp.tau_tdiff(:,:,:,band);
            
            ediff = interpol_2d_mex(tau_band, tau_ediff_band); % 8x1
            tdiff = interpol_3d_mex(tau_band, tau_tdiff_band); %9x8
            ediff_mix = fraction_band * ediff;
            tdiff_mix = tdiff * fraction_band';
            surf_lim(:,band) = (exp(-tau_band./smartp.mu) + tdiff_mix) * const.Config_albedo_thresh_land / pi * ...
                     (smartp.mu0 * exp(-tau_band/smartp.mu0) + ediff_mix);
        end  
       
    end

end
