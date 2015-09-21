function [resid,surf] = get_resid_eof2(regp,atm_path,const)

    resid = NaN * ones(const.Cam_Dim, const.Band_Dim);
    surf = NaN * ones(const.Cam_Dim,const.Band_Dim);

    cam_used = regp.channel_is_used(:,const.Band_Green);
    num_cam_used = sum(cam_used);
    diff = regp.mean_equ_ref(cam_used,:) - atm_path(cam_used, :);
    eof = reshape(regp.eof(1:num_cam_used, 1:num_cam_used),num_cam_used,num_cam_used);

    exp_coef = eof' * diff;
    idx = 1:num_cam_used <= regp.max_usable_eof;
    surf(cam_used, :) = eof(:,idx) * exp_coef(idx,:);
    resid(cam_used, :) = diff - surf(cam_used,:);
    
    resid = resid(:);
    surf = surf(:);

end