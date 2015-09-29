function [atm_path,surf,resid] = par_update_resid_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,reg,smart,tau,theta,ExtCroSect,CompSSA,const,kf,add_limit);
    
    sub_block_size = sub_block_end-sub_block_start+1;
    atm_path = NaN * ones(const.NChannel,sub_block_size);
    surf = NaN * ones(const.NChannel,sub_block_size);
    resid = NaN * ones(const.NChannel,sub_block_size);
    
    cnt = 1;
    for p =  sub_block_start : min(sub_block_end,num_pixel)

        xp = x(p);
        yp = y(p); 
        thetap = theta(:,p); 
        taup = tau(p);
        
        [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf,add_limit);
        [atm_path(:,cnt),surf(:,cnt),resid(:,cnt)] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
        
        cnt = cnt + 1;

    end
    
end


