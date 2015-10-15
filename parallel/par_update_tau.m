function [new_tau,new_resid] = par_update_tau(old_tau,theta,old_resid,kappa,sigmasq,delta, ...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, r, par, core,add_limit, const)
    
    if par == true
        
        num_pixel = reg.num_reg_used;
        size_sub_block = ceil(num_pixel/core);
        
        new_tau = NaN*ones(size_sub_block,core);
        new_resid = NaN*ones(const.NChannel,size_sub_block,core);
        
        parfor sub_block = 1:core
            
            sub_block_start = (sub_block-1)*size_sub_block+1;
            sub_block_end = sub_block*size_sub_block;
            
            [new_tau(:,sub_block),new_resid(:,:,sub_block)] = ...
            par_update_tau_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,old_resid,theta,kappa,sigmasq,delta,Method,ExtCroSect,CompSSA,i,j,...
            old_tau,reg,smart,r,add_limit,const);

        end
        
        new_tau = reshape(new_tau,size_sub_block*core,1);
        new_tau = new_tau(1:num_pixel);
        new_resid = reshape(new_resid,const.NChannel,size_sub_block*core);
        new_resid = new_resid(:,1:num_pixel);
        
    else
        
        for p = 1:reg.num_reg_used

            [old_residp,old_taup,thetap,old_tau_neighbor,regp,smartp] = par_preprocess_tau(x,y,p,old_resid,theta,old_tau,i,j,reg,smart,r,add_limit,const);

            [old_tau(p),old_resid(:,p)] = par_update_tau_pixel(old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
            regp,smartp,ExtCroSect,CompSSA,r,add_limit,const);        

        end
        
        new_tau = old_tau;
        new_resid = old_resid;
    end
        
    
end