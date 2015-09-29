function [atm_path,surf,resid] = par_update_resid(tau,theta, x, y, smart, reg, ExtCroSect, CompSSA, kf, par,core,add_limit, const)
   
    if par == true
        
        num_pixel = reg.num_reg_used;
        size_sub_block = ceil(num_pixel/core);
        
        atm_path = NaN*ones(const.NChannel,size_sub_block,core);
        surf = NaN*ones(const.NChannel,size_sub_block,core);
        resid = NaN*ones(const.NChannel,size_sub_block,core);
        
        parfor sub_block = 1:core
            
            sub_block_start = (sub_block-1)*size_sub_block+1;
            sub_block_end = sub_block*size_sub_block;
            
            [atm_path(:,:,sub_block),surf(:,:,sub_block),resid(:,:,sub_block)] = ...
            par_update_resid_sub_block(sub_block_start,sub_block_end,num_pixel,x,y,reg,smart,tau,theta,ExtCroSect,CompSSA,const,kf,add_limit);

        end
        
        atm_path = reshape(atm_path,const.NChannel,size_sub_block*core);
        atm_path = atm_path(:,1:num_pixel);
        surf = reshape(surf,const.NChannel,size_sub_block*core);
        surf = surf(:,1:num_pixel);
        resid = reshape(resid,const.NChannel,size_sub_block*core);
        resid = resid(:,1:num_pixel);
        
    else
        
        atm_path = NaN*ones(const.NChannel,reg.num_reg_used);
        surf = NaN*ones(const.NChannel,reg.num_reg_used);
        resid = NaN*ones(const.NChannel,reg.num_reg_used);

        for p = 1:reg.num_reg_used
            
            xp = x(p);
            yp = y(p); 
            thetap = theta(:,p); 
            taup = tau(p);

            [regp,smartp] = extract_pixel(xp,yp,reg,smart,const,kf,add_limit);

            [atm_path(:,p),surf(:,p),resid(:,p)] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
        end
    end

end