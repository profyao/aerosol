function show(sample,reg,band,cam,iter,cmap,const)
    
    tau = sample.tau(:,iter);
    atm_path = reshape(sample.atm_path((band-1)*9+cam,:,iter),reg.num_reg_used,1);
    surf = reshape(sample.surf((band-1)*9+cam,:,iter),reg.num_reg_used,1);
    resid = reshape(sample.resid((band-1)*9+cam,:,iter),reg.num_reg_used,1);
    ref = reshape(reg.mean_equ_ref(:,:,band,cam),const.XDim_r,const.YDim_r);
    [x,y] = find(reg.reg_is_used);
    
    if max(tau)~=min(tau)
        subplot(511)
        plot_1d(tau,x,y,cmap,const,[min(tau),max(tau)]);
        title('AOD')
    else
        subplot(511)
        plot_1d(tau,x,y,cmap,const);
        title('AOD')
    end
        
    subplot(512)
    plot_1d(atm_path,x,y,cmap,const);
    title('atm_path','interpreter','none')
    subplot(513)
    plot_1d(surf,x,y,cmap,const);
    title('surf')
    subplot(514)
    plot_1d(resid,x,y,cmap,const);
    title('residual')
    subplot(515)
    plot_2d(ref,cmap,const,[min(ref(:)),max(ref(:))]);
    title('reflectance')

end