function [resid,surf] = get_resid_eof(channel_is_used,mean_equ_ref,eof,max_usable_eof,atm_path,const)

    Cam_Dim = const.Cam_Dim;
    Band_Dim = const.Band_Dim;

    resid = NaN * ones(Cam_Dim, Band_Dim);
    surf = NaN * ones(Cam_Dim,Band_Dim);

    for band =1:Band_Dim

        cam_used = channel_is_used(:,band);
        num_cam_used = sum(cam_used);
        diff = mean_equ_ref(cam_used,band) - atm_path(cam_used, band);
        eof_band = reshape(eof(band, 1:num_cam_used, 1:num_cam_used),num_cam_used,num_cam_used);

        exp_coef = diff' * eof_band;
        idx = 1:num_cam_used <= max_usable_eof(band);
        surf(cam_used, band) = eof_band * (idx.*exp_coef)';
        resid(cam_used, band) = diff - surf(cam_used,band);

    end
    
    resid = resid(:);
    surf = surf(:);

end