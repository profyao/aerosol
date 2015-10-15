function [new_tau,new_resid] = par_update_tau_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,old_resid,theta,kappa,sigmasq,delta,Method,ExtCroSect,CompSSA,i,j,old_tau,reg,smart,r,add_limit,const)
    
    sub_block_size = sub_block_end-sub_block_start+1;
    new_tau = NaN * ones(sub_block_size,1);
    new_resid = NaN * ones(const.NChannel,sub_block_size);
    
    cnt = 1;
    for p =  sub_block_start : min(sub_block_end,num_pixel)

        [old_residp,old_taup,thetap,old_tau_neighbor,regp,smartp] = par_preprocess_tau(x,y,p,old_resid,theta,old_tau,i,j,reg,smart,r,add_limit,const);

        [new_tau(cnt),new_resid(:,cnt)] = par_update_tau_pixel(old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
        regp,smartp,ExtCroSect,CompSSA,r,add_limit,const);
        
        cnt = cnt + 1;

    end
    
end

