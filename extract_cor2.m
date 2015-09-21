function [cor1,cor2,ratio,component] = extract_cor2(reg,smart,x,y,ExtCroSect,CompSSA,kf,const)
    
    Model_ComponentDim = 21;
    theta_grid = eye(Model_ComponentDim);
    atm_path = NaN*ones(const.Cam_Dim,const.Band_Dim,const.Model_OpticalDepthLen,Model_ComponentDim);
    
    for i = 1:const.Model_OpticalDepthLen
        for j = 1:Model_ComponentDim
            tau = const.Model_OpticalDepthGrid(i);
            theta = theta_grid(:,j);
            atm_path(:,:,i,j) = get_model(tau,theta,x,y,ExtCroSect,CompSSA,smart,const);
        end
    end
    
    totdim = const.Model_OpticalDepthLen*Model_ComponentDim;
    atm_path = reshape(atm_path,const.Cam_Dim, const.Band_Dim,totdim);
    L = reshape(reg.mean_equ_ref(x,y,:,:),const.Band_Dim,const.Cam_Dim)';
    if kf == true
        top_eof = reshape(reg.eof(x,y,:,1),const.Cam_Dim,1);
    else
        top_eof = reshape(reg.eof(x,y,:,:,1),const.Band_Dim,const.Cam_Dim)';
    end
    
    cor1p = NaN*ones(const.Band_Dim,totdim);
    cor2 = NaN*ones(1,const.Band_Dim);
    ratiop = NaN*ones(const.Band_Dim,totdim);
    
    for band = 1:const.Band_Dim
        cor1p(band,:) = arrayfun(@(x) corr_nan(atm_path(:,band,x),L(:,band)), 1:totdim);
        if kf == true
            cor2(band) = corr_nan(top_eof,L(:,band));
        else
            cor2(band) = corr_nan(top_eof(:,band),L(:,band));
        ratiop(band,:) = arrayfun(@(x) norm(atm_path(:,band,x))/norm(L(:,band)), 1:totdim);
        end
    end

    cor1p = reshape(cor1p,[const.Band_Dim,const.Model_OpticalDepthLen,Model_ComponentDim]);
    ratiop = reshape(ratiop,[const.Band_Dim,const.Model_OpticalDepthLen,Model_ComponentDim]);
    
    cor1 = NaN*ones(1,const.Band_Dim);
    ratio = NaN*ones(1,const.Band_Dim);
    component = NaN*ones(1,const.Band_Dim);
    
    for band = 1:const.Band_Dim
       [tmp,id1] = max(reshape(cor1p(band,:,:),const.Model_OpticalDepthLen,Model_ComponentDim));
       [cor1(band),id2] = max(tmp);
       component(band) = id2;
       ratio(band) = ratiop(band,id1(id2),id2); 
    end
    
end