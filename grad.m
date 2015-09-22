function [g_v,flag] = grad(taup,tau_neighbor,thetap,var_str,sigmasq,residp,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit,varargin)

    switch var_str
        
        case 'tau'

            if isinf(residp(1))              
                g_v = 0;
                flag = -1; 
                
            else

                g_v = 0;
                flag = NaN;
                cnt = 1;
                max_iter = 10;
                delta_tau = 0.01;
                n_neighbor = length(tau_neighbor);
                kappa = varargin{1};

                while true

                    g0 = g_v;

                    taul = max(taup - delta_tau,0);
                    taur = min(taup + delta_tau,3);
                    [~,~,residl] = get_resid(taul,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
                    [~,~,residr] = get_resid(taur,thetap,regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);

                    if isinf(residl(1)) || isinf(residr(1))
                        g_v = 0;
                        flag = 0;
                        break
                    elseif n_neighbor > 0 
                          smoothl = kappa * sum(taul - tau_neighbor).^2;
                          smoothr = kappa * sum(taur - tau_neighbor).^2;
                    else
                      smoothl=0;smoothr=0;
                    end

                    chisql = nansum(residl(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
                    chisqr = nansum(residr(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));

                    g_v = (chisql-chisqr+smoothl-smoothr)/(taul-taur);
                    
                    if abs(g_v-g0)<abs(g0)*1e-3
                        break
                    end
                    
                    if cnt > max_iter
                        g_v = 0;
                        flag = 0;
                        fprintf('cannot converge to gradient within %d iterations! delta %s: %e!\n',max_iter,var_str,delta_tau)
                        break
                    end
                    
                    delta_tau = 0.5 * delta_tau;
                    cnt = cnt + 1;
                end
                
                if isnan(flag)
                    flag = 1;
                end

            end

        case 'theta'  
            
            if isinf(residp(1))
                g_v = zeros(const.Component_Num,1);
                flag = -1 * ones(const.Component_Num,1);               
            else
                
                g_v = NaN*ones(const.Component_Num,1);
                flag = NaN*ones(const.Component_Num,1);
           
                g = 0; 
                cnt = 1;
                max_iter = 10;
              
                for c = 1:const.Component_Num
                    
                    delta_theta = 0.001;
                    delta_theta_v = zeros(const.Component_Num,1);
                    delta_theta_v(c) = delta_theta;
                    
                    while true                   
                        
                        g0 = g;

                        thetal = thetap - delta_theta_v; thetal(thetal<0)=0;
                        thetar = thetap + delta_theta_v;
                        [~,~,residl] = get_resid(taup,thetal/sum(thetal),regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);
                        [~,~,residr] = get_resid(taup,thetar/sum(thetar),regp,smartp,ExtCroSect,CompSSA,const,kf,add_limit);

                        if isinf(residl(1)) || isinf(residr(1))
                            g_v(c) = 0;
                            flag(c) = 0;
                            break 
                            
                        else

                            chisql = nansum(residl(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));       
                            chisqr = nansum(residr(const.Channel_Used).^2 ./ sigmasq(const.Channel_Used));

                            g = (chisql-chisqr)/(thetal(c)-thetar(c));
                            
                            if abs(g-g0)<abs(g0)*1e-3
                                break
                            end   
                            
                            if cnt > max_iter
                                g_v(c) = 0;
                                flag(c) = 0;
                                fprintf('cannot converge to gradient within %d iterations! delta %s: %e!\n',max_iter,var_str,delta_theta)
                                break
                            end

                            delta_theta = 0.5 * delta_theta;
                            cnt = cnt + 1;

                        end
                    end
                    
                    if isnan(flag(c))
                        g_v(c) = g;
                        flag(c) = 1;
                    end

                end
            end

    end
            
end