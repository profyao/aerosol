function [cor1,cor2,ratio] = extract_cor(reg,smart,scatter_type,x,y,band_ind,const)

    x0 = ceil(x/const.RegScale);
    y0 = ceil(y/const.RegScale);
    
    if strcmp(scatter_type,'SS')
        rho = permute(reshape(smart.ss(:,x0,y0,:,band_ind,:),const.Model_OpticalDepthLen,const.Model_ComponentDim,const.Cam_Dim)...
            ,[3,1,2]);
    elseif strcmp(scatter_type,'MS')
        rho = permute(reshape(smart.ms(:,x0,y0,:,band_ind,:),const.Model_OpticalDepthLen,const.Model_ComponentDim,const.Cam_Dim)...
            ,[3,1,2]);
    else
        error('need to specify scattering type!\n')
    end
    
    totdim = const.Model_OpticalDepthLen*const.Model_ComponentDim;

    rho = reshape(rho,const.Cam_Dim,totdim);
    L = reshape(reg.mean_equ_ref(x,y,band_ind,:),const.Cam_Dim,1);
    top_eof = reshape(reg.eof(x,y,:,1),const.Cam_Dim,1);

    cor1 = arrayfun(@(x) corr(rho(:,x),L), 1:totdim);
    cor2 = nancov(top_eof,L)/(nanstd(top_eof)*nanstd(L));
    ratio = arrayfun(@(x) norm(rho(:,x))/norm(L), 1:totdim);
    
    cor1 = reshape(cor1,[const.Model_OpticalDepthLen,const.Model_ComponentDim]);
    ratio = reshape(ratio,[const.Model_OpticalDepthLen,const.Model_ComponentDim]);

    
end