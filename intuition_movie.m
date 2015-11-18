function intuition_movie(Date,Path,Orbit,Block,r,Method)

    [reg, sample] = load_cache(Date,Path,Orbit,Block,r,'reg','sample',Method);
    [x,y] = find(reg.reg_is_used);
    
    %for p = 1:reg.num_reg_used
        
    %end
    
    p = 100;
    atm_path = sample.atm_path(:,p);
    surf = sample.surf(:,p);
    resid = sample.resid(:,p);
    sigma = sqrt(sample.sigmasq);
    
    theta = sample.theta(:,p)
    
    obs = squeeze(reg.mean_equ_ref(x(p),y(p),:,:))';
    obs = obs(:);
    
    figure
    plot([obs, atm_path, surf, resid],'-','LineWidth',2)
    legend({'obs','atm\_path','surf','resid'},'Location','northwest')
    
    set(gca,'FontSize',18)

end