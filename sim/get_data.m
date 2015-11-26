function reg_sim = get_data(r, XDim_r, YDim_r, const, reg, smart, x, y, tau, theta, ExtCroSect,CompSSA)

    surf = zeros(const.Cam_Dim,const.Band_Dim);
    %surf = kron([0.05, 0.05, 0.05, 0.5],cosd(degree));
    noise = zeros(const.Cam_Dim,const.Band_Dim);

    reg_sim = NaN*ones(XDim_r,YDim_r,const.Band_Dim,const.Cam_Dim);
    for p = 1:reg.num_reg_used
        xp = x(p);
        yp = y(p);
        taup = tau(p);
        thetap = theta(:,p);
        [~,smartp] = extract_pixel(xp,yp,reg,smart,const,r,0);
        [atm_path,~] = get_model(taup,thetap,ExtCroSect,CompSSA,smartp,const,0);
        equ_refp = atm_path + surf + noise;
        reg_sim(xp,yp,:,:) = equ_refp';
    end
    
    figure,plot_2d(reg_sim(:,:,2,1),jet),colorbar
end