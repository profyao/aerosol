function [new_taup,new_residp] = par_update_tau_pixel(old_residp,old_taup,old_tau_neighbor,thetap,kappa,sigmasq,delta,Method,...
    regp,smartp,ExtCroSect,CompSSA, kf, add_limit, const)

        n_neighbor = length(old_tau_neighbor);
        if n_neighbor > 0
            mean_neighbor = mean(old_tau_neighbor);
        end
        
        if strcmp(Method,'CD-random') || strcmp(Method,'MCMC') || strcmp(Method,'CD-random-noprior')
   
            if isinf(old_residp(1))
                taup = old_taup * 0.8;
                smooth0 = 0;
                smooth1 = 0;
            else

                if n_neighbor > 0                   
                    mu = 0.5 * (mean_neighbor + old_taup);
                    taup = mu + delta * randn(1);
                    taup(taup<0)=0;
                    taup(taup>3)=3;
                    smooth0 = kappa * sum(old_taup - old_tau_neighbor).^2;
                    smooth1 = kappa * sum(taup-old_tau_neighbor).^2;
                else
                    mu = old_taup;
                    taup = mu + delta * randn(1);
                    taup(taup<0)=0;
                    taup(taup>3)=3;
                    smooth0=0;smooth1=0;
                end     

            end
            
            [~,~,residp] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);            

            new_chisq = nansum(residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
            chisq = nansum(old_residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));

            if isinf(old_residp(1)) && isinf(residp(1))
                new_taup = taup;
                new_residp = residp;
                return
            end

            switch Method

                case {'CD-random','CD-random-noprior'}

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
            
        elseif strcmp(Method,'CD')
            
            [g,flag] = grad(old_taup,old_tau_neighbor,thetap,'tau',sigmasq,old_residp,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit,kappa);
            
            if flag == -1
                
                new_taup = old_taup * 0.8;
                [~,~,new_residp] = get_resid(new_taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
                           
            elseif flag == 0
                
                new_taup = old_taup;
                new_residp = old_residp;
                
            elseif flag == 1
                
                [new_residp,new_taup] = back_track(g,old_taup,old_tau_neighbor,thetap,'tau',sigmasq,old_residp,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit,kappa);                            
                
            else          
                error('flag is not assigned when calculating gradient!')
            end            
            
        else
            error('No Method is specified!\n')
        end
 
end