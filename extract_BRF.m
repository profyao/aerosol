function [theta0,theta,omega,phi,rho] = extract_BRF(aod_ind,srf_pres_ind,band_ind,component_ind,mu0_ind,scatter_type,const)
    
    if strcmp(scatter_type,'SS')
        file_smart = fullfile('products/MIANSMT/',const.MIANSMT_SS_filename);
    elseif strcmp(scatter_type,'MS')
        file_smart = fullfile('products/MIANSMT/',const.MIANSMT_MS_filename);
    else
        error('incorrect scattering type!\n')
    end

    EquivalentReflectanceMu0 = hdfread(file_smart,['/EquivalentReflectanceMu0_',num2str(mu0_ind)], 'Index', ...
            {[1  1  1  1  1  1],[1  1  1  1  1  1],[4  21  13  2  29  96]});
    
    EquivalentReflectanceMu0 = reshape(EquivalentReflectanceMu0(band_ind,component_ind,aod_ind,srf_pres_ind,:,:),29,96);
    
    mu0_grid = const.Model_mu0Grid;
    mu0 = mu0_grid(mu0_ind);
    theta0 = acosd(mu0);
    
    mu_grid = const.Model_muGrid;
    omega_grid = const.Model_ScatterAngleGrid;

    [mu_ind,omega_ind] = find(EquivalentReflectanceMu0~=32767);
    rho = EquivalentReflectanceMu0(EquivalentReflectanceMu0~=32767);
    rho = exp(double(rho)/1000);
    
    mu = mu_grid(mu_ind);
    theta = acosd(mu)';
    omega = omega_grid(omega_ind)';
    
    theta_grid = acosd(mu_grid);
    low_allowable_omega = 180 - (theta_grid + theta0);
    high_allowable_omega = 180 - abs(theta_grid-theta0);
    omega(omega==-1) = low_allowable_omega;
    omega(omega==181) = high_allowable_omega;
    
    omega = 180-omega;
    %phi = NaN* ones(length(theta),1);
   
    tmp1 = cosd(omega) - cosd(theta0)*cosd(theta);
    tmp2 = sind(theta0)*sind(theta);
    %idx = tmp1 == 0 & tmp2 ==0;

    cosd_phi = tmp1./tmp2;
    cosd_phi(cosd_phi>1) = 1;
    cosd_phi(cosd_phi<-1) = -1;

    phi = acosd(cosd_phi);
    
end