function [new_thetap,new_residp] = par_update_theta_pixel(old_residp,taup,old_thetap,old_theta_neighbor,sigmasq,alpha,Method,...
            regp,smartp,ExtCroSect,CompSSA,r,add_limit,const)
    
        n_neighbor = length(old_theta_neighbor);
        
        if strcmp(Method,'CD-random') || strcmp(Method,'CD-random-noprior') || strcmp(Method,'MCMC-G')
        
            if n_neighbor > 0
                mu = mean(old_theta_neighbor,2);
            else
                mu = old_thetap;
            end
            
            thetap = gamrnd(mu,1);
            %thetap(thetap<1e-6) = 1e-6;
            %thetap = gamrnd(0.05*ones(const.Component_Num,1),1);
            thetap = thetap / sum(thetap);

            [~,~,residp] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit);
            
            new_chisq = nansum(residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
            chisq = nansum(old_residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));
            
            if strcmp(Method,'MCMC-G')
                l1 = (alpha-1)'*log(thetap);
                l0 = (alpha-1)'*log(old_thetap);
                
                if isinf(old_residp(1)) && isinf(residp(1))
                    new_thetap = thetap;
                    new_residp = residp;  
                elseif chisq - l0 > new_chisq - l1
                    new_thetap = thetap;
                    new_residp = residp;
                else
                    new_thetap = old_thetap;
                    new_residp = old_residp;  
                end
            else
                
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
            end
                                                       
        elseif strcmp(Method,'MCMC')
        
            thetap = gamrnd(alpha, 1); % Dirichlet independent proposal
            thetap = thetap/sum(thetap);
        
            [~,~,residp] = get_resid(taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit);
            
            if isinf(old_residp(1)) && isinf(residp(1))
                new_thetap = old_thetap;
                new_residp = old_residp;
            else
                
                new_chisq = nansum(residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
                chisq = nansum(old_residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used)); 
                
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
            
        elseif strcmp(Method,'CD')
            
            [g,~] = grad(taup,old_theta_neighbor,old_thetap,'theta',sigmasq,old_residp,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit);
                            
            [new_residp,new_thetap] = back_track(g,taup,old_theta_neighbor,old_thetap,'theta',sigmasq,old_residp,regp,smartp,ExtCroSect,CompSSA,const,r,add_limit);
            
        else
            
            error('No Method is specified!\n')
            
        end

end