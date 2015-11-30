function chisq_grid = gen_tt_grid(r,theta_misr_grid,sigmasq,xp,yp,ExtCroSect,CompSSA,smart,reg,const)
    
    chisq_grid = zeros(13,74);
    for ii = 1:13
        for jj = 1:74
            tau_tmp = const.Model_OpticalDepthGrid(ii);
            theta_tmp = theta_misr_grid(:,jj);
            chisq_grid(ii,jj) = 0.5 * get_chisq(r,tau_tmp,theta_tmp,sigmasq,xp,yp,const,ExtCroSect,CompSSA,smart,reg);
        end
    end

end