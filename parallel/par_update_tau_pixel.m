function [new_taup,new_residp] = par_update_tau_pixel(xp,yp,old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
    channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,ExtCroSect,CompSSA,smart,const)


        n_neighbor = length(old_tau_neighbor);
        if n_neighbor > 0
            mean_neighbor = mean(old_tau_neighbor);
        end
        
        if strcmp(Method,'CDSS') || strcmp(Method,'MCMC')
   
            if isinf(old_residp(1))
                taup = old_taup * 0.8;
                smooth0 = 0;
                smooth1 = 0;
            else

                if n_neighbor > 0                   
                    mu = 0.5 * (mean_neighbor + old_taup);
                    taup = mu + delta * randn(1);
                    taup(taup<=1e-3)=1e-3;
                    taup(taup>=3)=3;
                    smooth0 = kappa * sum(old_taup - old_tau_neighbor).^2;
                    smooth1 = kappa * sum(taup-old_tau_neighbor).^2;
                else
                    mu = old_taup;
                    taup = mu + delta * randn(1);
                    taup(taup<=1e-3)=1e-3;
                    taup(taup>=3)=3;
                    smooth0=0;smooth1=0;
                end     

            end

            [~,~,residp] = get_resid(taup,thetap,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);

            new_chisq = nansum(residp.^2 ./ sigmasq);       
            chisq = nansum(old_residp.^2 ./ sigmasq);

            if isinf(old_residp(1)) && isinf(residp(1))
                new_taup = taup;
                new_residp = residp;
                return
            end

            switch Method

                case 'CDSS'

                    if chisq + smooth0 > new_chisq + smooth1
                        new_taup = taup;
                        new_residp = residp;
                    else
                        new_taup = old_taup;
                        new_residp = old_residp;
                    end

                case 'MCMC'

                    if n_neighbor > 0
                        term1 = 0.5 * (sum((taup - old_tau_neighbor).^2) - sum((old_taup - old_tau_neighbor).^2));
                        term2 = 0.5 * (new_chisq - chisq);
                        term3 = 0.5 * n_neighbor / delta * ( (taup-0.5*(old_taup + mean_neighbor))^2 - (old_taup - 0.5*(taup+mean_neighbor))^2 );
                        A = exp( - term1 - term2 + term3);
                    else
                        A = exp( - 0.5* (new_chisq - chisq) );
                    end

                    u = rand(1);

                    if u < A 
                        new_taup = taup;
                        new_residp = residp;
                    else
                        new_taup = old_taup;
                        new_residp = old_residp;
                    end    
            end
            
        elseif strcmp(Method,'MAP')
            
            [g,flag] = grad(old_taup,old_tau_neighbor,thetap,'tau',sigmasq,old_residp,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const,kappa);
            
            if flag == -1
                
                new_taup = old_taup * 0.8;
                [~,~,new_residp] = get_resid(new_taup,thetap,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const);
                           
            elseif flag == 0
                
                new_taup = old_taup;
                new_residp = old_residp;
                
            elseif flag == 1
                
                [new_residp,new_taup] = back_track(g,old_taup,old_tau_neighbor,thetap,'tau',sigmasq,old_residp,xp,yp,channel_is_used,min_equ_ref,mean_equ_ref,eof,max_usable_eof,...
                smart,ExtCroSect,CompSSA,const,kappa);                            
                
            else          
                error('flag is not assigned when calculating gradient!')
            end            
            
        else
            error('No Method is specified!\n')
        end
 
end