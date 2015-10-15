function [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,r,add_limit)

    regp.channel_is_used = reshape(reg.channel_is_used(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';
    regp.min_equ_ref = reshape(reg.min_equ_ref(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';
    regp.mean_equ_ref = reshape(reg.mean_equ_ref(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';
    
    RegScale = const.r17600/r;

    if r > 1100
        regp.eof = reshape(reg.eof(xp,yp,:,:,:),const.Band_Dim,const.Cam_Dim,const.Cam_Dim);
        regp.max_usable_eof = reshape(reg.max_usable_eof(xp,yp,:),const.Band_Dim,1);
    else
        regp.eof = reshape(reg.eof(xp,yp,:,:),const.Cam_Dim,const.Cam_Dim);
        regp.max_usable_eof = reg.max_usable_eof(xp,yp);
    end

    smartp.tau_cam_ss = permute(reshape(smart.ss(:, ceil(xp/RegScale),ceil(yp/RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);
    smartp.tau_cam_ms = permute(reshape(smart.ms(:, ceil(xp/RegScale),ceil(yp/RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);
    
    if add_limit == true
        smartp.ss_ediff = reshape(smart.ss_ediff(:, ceil(xp/RegScale), ceil(yp/RegScale), const.Component_Particle, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim);
        smartp.ms_ediff = reshape(smart.ms_ediff(:, ceil(xp/RegScale), ceil(yp/RegScale), const.Component_Particle, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim);

        smartp.ss_tdiff = permute(reshape(smart.ss_tdiff(:, ceil(xp/RegScale), ceil(yp/RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);
        smartp.ms_tdiff = permute(reshape(smart.ms_tdiff(:, ceil(xp/RegScale), ceil(yp/RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);

        smartp.mu = reshape(smart.mu(ceil(xp/RegScale),ceil(yp/RegScale),:),const.Cam_Dim,1);
        smartp.mu0 = smart.mu0(ceil(xp/RegScale), ceil(yp/RegScale));
    end


end