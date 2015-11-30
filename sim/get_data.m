function [reg_sim,noise] = get_data(r, XDim_r, YDim_r, const, reg, smart, x, y, tau, theta, ExtCroSect,CompSSA, surf, ratio)

    reg_sim = NaN*ones(XDim_r,YDim_r,const.Band_Dim,const.Cam_Dim);
    noise = NaN*ones(XDim_r,YDim_r,const.Band_Dim,const.Cam_Dim);
    
    for p = 1:reg.num_reg_used
        xp = x(p);
        yp = y(p);
        taup = tau(p);
        thetap = theta(:,p);
        [~,smartp] = extract_pixel(xp,yp,reg,smart,const,r,0);
        [atm_path,~] = get_model(taup,thetap,ExtCroSect,CompSSA,smartp,const,0);
        avg = mean(atm_path);
        sd = ratio*avg;
        noisep = repmat(sd,const.Cam_Dim,1) .* randn(const.Cam_Dim,const.Band_Dim);
        equ_refp = atm_path + surf + noisep;
        reg_sim(xp,yp,:,:) = equ_refp';
        noise(xp,yp,:,:) = noisep';
    end
    
    %figure,plot_2d(reg_sim(:,:,2,1),jet),colorbar
end