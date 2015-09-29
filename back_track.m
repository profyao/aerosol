function [new_residp,new_var] = back_track(g,taup,tau_neighbor,thetap,var_str,sigmasq,residp,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit,varargin)
    
      
    switch var_str
        
        case 'tau'
            
            if ~any(g)
                new_var = taup;
                new_residp = residp;                
            else
                                
                kappa = varargin{1};
                cnt = 1;
                max_iter = 10;
                lambda = 1e-4;
                n_neighbor = length(tau_neighbor);

                while true

                    new_taup = max(taup - lambda*g,1e-3);
                    new_taup = min(new_taup,3);
                    
                    [~,~,new_residp] = get_resid(new_taup,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);

                    if n_neighbor > 0 
                        new_smooth = kappa * sum(new_taup - tau_neighbor).^2;
                        smooth = kappa * sum(taup - tau_neighbor).^2;
                    else
                        new_smooth=0;smooth=0;
                    end

                    new_chisq = nansum(new_residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
                    chisq = nansum(residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));

                    if new_chisq+new_smooth < chisq+smooth
                        new_var = new_taup;
                        break
                    end

                    if cnt > max_iter
                        new_var = taup;
                        new_residp = residp;
                        fprintf('cannot achieve descent within %d iterations! %s: %e, last value: %e, curret value: %e!\n',max_iter,var_str,lambda,chisq+smooth,new_chisq+new_smooth);
                        break
                    end

                    lambda = 0.5 * lambda;
                    cnt = cnt + 1;

                end
            end
                                  
        case 'theta'
            
            if ~any(g)
                new_var = thetap;
                new_residp = residp;
            else
            
                cnt = 1;
                max_iter = 10;
                lambda = 1e-4;

                while true

                    new_thetap = thetap - lambda*g;
                    new_thetap(new_thetap<0) = 0;
                    new_thetap = new_thetap/sum(new_thetap);

                    [~,~,new_residp] = get_resid(taup,new_thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);


                    new_chisq = nansum(new_residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
                    chisq = nansum(residp(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));

                    if new_chisq < chisq
                        new_var = new_thetap;
                        break
                    end

                    if cnt > max_iter
                        new_var = thetap;
                        new_residp = residp;
                        fprintf('cannot achieve descent within %d iterations! %s: %e, last value: %e, curret value: %e!\n',max_iter,var_str,lambda,chisq,new_chisq);
                        break
                    end

                    lambda = 0.5 * lambda;
                    cnt = cnt + 1;

                end
            end
                      
    end

end