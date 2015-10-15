function [cor1,cor2,cor3,ratio] = extract_cor(reg,sample,xid,yid,const)
    
    valid = find(reg.reg_is_used);
    [x,y] = ind2sub([const.XDim_r,const.YDim_r],valid);
    id = ismember([x,y],[xid,yid],'rows');
    tmp = reshape(reg.mean_equ_ref,[const.XDim_r*const.YDim_r,const.Band_Dim,const.Cam_Dim]);
    tmp = reshape(tmp(valid,1,:),reg.num_reg_used,const.Cam_Dim)';
    L = tmp(:,id); 
    atm_path = squeeze(sample.atm_path(1:const.Cam_Dim,id,end));
    surf = squeeze(sample.surf(1:const.Cam_Dim,id,end));
    resid = squeeze(sample.resid(1:const.Cam_Dim,id,end));
    num = sum(id);
    cor1 = arrayfun(@(x) corr_nan(L(:,x),atm_path(:,x)), 1:num);
    cor2 = arrayfun(@(x) corr_nan(atm_path(:,x),surf(:,x)), 1:num);
    cor3 = arrayfun(@(x) corr_nan(L(:,x),resid(:,x)), 1:num);
    ratio = arrayfun(@(x) norm_nan(atm_path(:,x))/norm_nan(L(:,x)), 1:num);
    
end