function [new_thetap,new_residp] = par_update_theta_pixel(xp,yp,old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
            channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,smart,ExtCroSect,CompSSA,const)
    
        n_neighbor = length(old_theta_neighbor);
        
        if strcmp(Method,'CDSS')
        
            if n_neighbor > 0
                mu = mean(old_theta_neighbor,2);
            else
                mu = old_thetap;
            end

            thetap = gamrnd(mu,1);
            thetap = thetap / sum(thetap);

            [~,~,residp] = get_resid(taup,thetap,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);
            
            new_chisq = nansum(residp.^2 ./ sigmasq);       
            chisq = nansum(old_residp.^2 ./ sigmasq);
            
            if isinf(old_residp(1)) && isinf(residp(1))
                new_thetap = thetap;
                new_residp = residp;  
            elseif chisq > new_chisq
                new_thetap = thetap;
                new_residp = residp;
            else
                new_thetap = old_thetap;
                new_residp = old_residp;
            end
                    
        elseif strcmp(Method,'MCMC')
        
            thetap = gamrnd(alpha, 1); % Dirichlet independent proposal
            thetap = thetap/sum(thetap);
        
            [~,~,residp] = get_resid(taup,thetap,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);
            
            if isinf(old_residp(1)) && isinf(residp(1))
                new_thetap = old_thetap;
                new_residp = old_residp;
            else
                
                new_chisq = nansum(residp.^2 ./ sigmasq);       
                chisq = nansum(old_residp.^2 ./ sigmasq);    

                A = exp(- 0.5 * (new_chisq - chisq));
                u = rand(1);

                if u < A
                    new_thetap = thetap;
                    new_residp = residp;
                else
                    new_thetap = old_thetap;
                    new_residp = old_residp;
                end
            end
            
        elseif strcmp(Method,'MAP')
            
            [g,~] = grad(taup,old_theta_neighbor,old_thetap,'theta',sigmasq,old_residp,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);
                            
            [new_residp,new_thetap] = back_track(g,taup,old_theta_neighbor,old_thetap,'theta',sigmasq,old_residp,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);
            
        else
            
            error('No Method is specified!\n')
            
        end

end