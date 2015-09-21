function [resid,surf] = get_resid_eof(regp,atm_path,const)

    resid = NaN * ones(const.Cam_Dim, const.Band_Dim);
    surf = NaN * ones(const.Cam_Dim,const.Band_Dim);

    for band =1:const.Band_Dim

        cam_used = regp.channel_is_used(:,band);
        num_cam_used = sum(cam_used);
        diff = regp.mean_equ_ref(cam_used,band) - atm_path(cam_used, band);
        eof_band = reshape(regp.eof(band, 1:num_cam_used, 1:num_cam_used),num_cam_used,num_cam_used);

        exp_coef = eof_band' * diff;
        idx = 1:num_cam_used <= regp.max_usable_eof(band);
        surf(cam_used, band) = eof_band(:,idx) * exp_coef(idx);
        resid(cam_used, band) = diff - surf(cam_used,band);

    end
    
    resid = resid(:);
    surf = surf(:);

end