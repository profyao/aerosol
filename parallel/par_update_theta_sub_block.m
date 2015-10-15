function [new_theta,new_resid] = par_update_theta_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,old_resid,old_theta,sigmasq,alpha,Method,ExtCroSect,CompSSA,i,j,tau,reg,smart,r,add_limit,const)
        
    sub_block_size = sub_block_end-sub_block_start+1;
    new_theta = NaN * ones(const.Component_Num,sub_block_size);
    new_resid = NaN * ones(const.NChannel,sub_block_size);
            
    cnt = 1;
    for p =  sub_block_start : min(sub_block_end,num_pixel)

        [old_residp,old_thetap,taup,old_theta_neighbor,regp,smartp] = par_preprocess_theta(x,y,p,old_resid,old_theta,tau,i,j,reg,smart,r,add_limit,const);

        [new_theta(:,cnt),new_resid(:,cnt)] = par_update_theta_pixel(old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
        regp,smartp,ExtCroSect,CompSSA,r,add_limit,const); 
    
        cnt = cnt + 1;

    end
    
end

