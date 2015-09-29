function [new_theta, new_resid] = par_update_theta(old_theta,tau,old_resid,sigmasq,alpha,...
    i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, kf, par, core, add_limit,const)
    
    if par == true
        
        num_pixel = reg.num_reg_used;
        size_sub_block = ceil(num_pixel/core);
        
        new_theta = NaN*ones(const.Component_Num,size_sub_block,core);
        new_resid = NaN*ones(const.NChannel,size_sub_block,core);
        
        parfor sub_block = 1:core
            
            sub_block_start = (sub_block-1)*size_sub_block+1;
            sub_block_end = sub_block*size_sub_block;
            
            [new_theta(:,:,sub_block),new_resid(:,:,sub_block)] = ...
            par_update_theta_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,old_resid,old_theta,sigmasq,alpha,Method,ExtCroSect,CompSSA,i,j,...
            tau,reg,smart,kf,add_limit,const);

        end
        
        new_theta = reshape(new_theta,const.Component_Num,size_sub_block*core);
        new_theta = new_theta(:,1:num_pixel);
        new_resid = reshape(new_resid,const.NChannel,size_sub_block*core);
        new_resid = new_resid(:,1:num_pixel);
        
    else
    
        for p = 1:reg.num_reg_used

            [old_residp,old_thetap,taup,old_theta_neighbor,regp,smartp] = par_preprocess_theta(x,y,p,old_resid,old_theta,tau,i,j,reg,smart,kf,add_limit,const);

            [old_theta(:,p),old_resid(:,p)] = par_update_theta_pixel(old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
                regp,smartp,ExtCroSect,CompSSA,kf,add_limit,const);   
        end
        
        new_theta = old_theta;
        new_resid = old_resid;
             
    end

end
