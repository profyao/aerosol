function show(r,current,reg,band,cam,cmap,const)
    
    XDim_r = const.XDim_r4400 * const.r4400/r;
    YDim_r = const.YDim_r4400 * const.r4400/r;
    
    tau = current.tau;
    atm_path = reshape(current.atm_path((band-1)*9+cam,:),reg.num_reg_used,1);
    surf = reshape(current.surf((band-1)*9+cam,:),reg.num_reg_used,1);
    resid = reshape(current.resid((band-1)*9+cam,:),reg.num_reg_used,1);
    ref = reshape(reg.mean_equ_ref(:,:,band,cam),XDim_r,YDim_r);
    [x,y] = find(reg.reg_is_used);    
    
    if max(tau)~=min(tau)
        subplot(511)
        plot_1d(r,tau,x,y,cmap,const,[0.2,0.35]);
        title('AOD')
    else
        subplot(511)
        plot_1d(r,tau,x,y,cmap,const);
        title('AOD')
    end
        
    subplot(512)
    plot_1d(r,atm_path,x,y,cmap,const);
    title('atm_path','interpreter','none')
    subplot(513)
    plot_1d(r, surf,x,y,cmap,const);
    title('surf')
    subplot(514)
    plot_1d(r, resid,x,y,cmap,const);
    title('residual')
    subplot(515)
    plot_2d(ref,cmap,const,[min(ref(:)),max(ref(:))]);
    title('reflectance')

end