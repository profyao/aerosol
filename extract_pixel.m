function [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf)

    regp.channel_is_used = reshape(reg.channel_is_used(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';
    regp.min_equ_ref = reshape(reg.min_equ_ref(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';
    regp.mean_equ_ref = reshape(reg.mean_equ_ref(xp,yp,:,:),const.Band_Dim,const.Cam_Dim)';

    if kf == false
        regp.eof = reshape(reg.eof(xp,yp,:,:,:),const.Band_Dim,const.Cam_Dim,const.Cam_Dim);
        regp.max_usable_eof = reshape(reg.max_usable_eof(xp,yp,:),const.Band_Dim,1);
    else
        regp.eof = reshape(reg.eof_allband(xp,yp,:,:),Cam_Dim,Cam_Dim);
        regp.max_usable_eof = reg.max_usable_eof_allband(xp,yp);
    end

    smartp.tau_cam_ss = permute(reshape(smart.ss(:, ceil(xp/const.RegScale),ceil(yp/const.RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);
    smartp.tau_cam_ms = permute(reshape(smart.ms(:, ceil(xp/const.RegScale),ceil(yp/const.RegScale), const.Component_Particle, :, :),const.Model_OpticalDepthLen,const.Component_Num,const.Band_Dim,const.Cam_Dim),[1,4,2,3]);
    smartp.tau_ediff = smart.ss_ediff(:, ceil(xp/const.RegScale), ceil(yp/const.RegScale), const.Component_Particle, :) + ...
            smart.ms_ediff(:, ceil(xp/const.RegScale), ceil(yp/const.RegScale), const.Component_Particle, :);
    smartp.tau_ediff = reshape(smartp.tau_ediff,const.Model.OpticalDepthLen,const.Component_Num,const.Band_Dim); 

    smartp.tau_tdiff = smart.ss_tdiff(:, ceil(xp/const.RegScale), ceil(yp/const.RegScale), const.Component_Particle, :, :) + ...
            smart.ms_tdiff(:, ceil(xp/const.RegScale), ceil(yp/const.RegScale), const.Component_Particle, :, :);
    smartp.tau_tdiff = permute(reshape(smartp.tau_tdiff,const.Model.OpticalDepthLen, const.Component_Num, const.Band_Dim, const.Cam_Dim),[1,4,2,3]);

    smartp.mu = reshape(smart.mu(ceil(xp/const.RegScale),ceil(yp/const.RegScale),:),const.Cam.Dim,1);
    smartp.mu0 = smart.mu0(ceil(xp/const.RegScale), ceil(yp/const.RegScale));


end