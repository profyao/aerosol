function [atm_path,surf_lim] = get_model(tau,theta,ExtCroSect,CompSSA,smartp,const,add_limit)
    
    atm_path = zeros(const.Cam_Dim,const.Band_Dim);
    surf_lim = zeros(const.Cam_Dim,const.Band_Dim);
    
    for band = 1:const.Band_Dim
        
        scale_factor = ExtCroSect(band, const.Component_Particle) ./ ExtCroSect(const.Band_Green, const.Component_Particle) * theta;
        fraction_band = ExtCroSect(band, const.Component_Particle)./ ExtCroSect(const.Band_Green, const.Component_Particle).* ...
            theta' / scale_factor;

        ssa_v = CompSSA(band, const.Component_Particle); %1x8
        ssa_mixture = ssa_v*fraction_band';

        tau_band = double(tau*scale_factor);
        ss_band = smartp.tau_cam_ss(:,:,:,band);
        ms_band = smartp.tau_cam_ms(:,:,:,band);
        rayleigh = reshape(ms_band(1,:,:),const.Cam_Dim,const.Component_Num);

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
                        
            ss_ediff_band = reshape(smartp.ss_ediff(:,:,band),const.Model_OpticalDepthLen,const.Component_Num);
            ms_ediff_band = reshape(smartp.ms_ediff(:,:,band),const.Model_OpticalDepthLen,const.Component_Num);
            ss_tdiff_band = reshape(smartp.ss_tdiff(:,:,:,band),const.Model_OpticalDepthLen,const.Cam_Dim,const.Component_Num);
            ms_tdiff_band = reshape(smartp.ms_tdiff(:,:,:,band),const.Model_OpticalDepthLen,const.Cam_Dim,const.Component_Num);
            
            ediff_band = ss_ediff_band + ms_ediff_band;
            tdiff_band = ss_tdiff_band + ms_tdiff_band;
            
            if tau_band == 0
                ediff = ediff_band(1,:)';
                tdiff = reshape(tdiff_band(1,:,:),const.Cam_Dim,const.Component_Num);
            else          
                ediff = interpol_2d_mex(tau_band, ediff_band); % 8x1
                tdiff = interpol_3d_mex(tau_band, tdiff_band); %9x8          
            end
            ediff_mix = fraction_band * ediff;
            tdiff_mix = tdiff * fraction_band';
            surf_lim(:,band) = (exp(-tau_band./smartp.mu) + tdiff_mix) * const.Config_albedo_thresh_land / pi * ...
                     (smartp.mu0 * exp(-tau_band/smartp.mu0) + ediff_mix);
        end        
        
    end

end