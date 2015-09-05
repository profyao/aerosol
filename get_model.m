function atm_path = get_model(tau,theta,x,y,ExtCroSect,CompSSA,smart,const)
    
    Cam_Dim = const.Cam_Dim;
    Band_Dim = const.Band_Dim;
    Band_Green = const.Band_Green;
    Component_Particle = const.Component_Particle;
    Component_Num = const.Component_Num;
    Model_OpticalDepthLen = const.Model_OpticalDepthLen;
    RegScale = const.RegScale;
    
    atm_path = zeros(Cam_Dim,Band_Dim);
    
    for band = 1:Band_Dim
        
        scale_factor = ExtCroSect(band, Component_Particle) ./ ExtCroSect(Band_Green, Component_Particle) * theta;
        fraction_band = ExtCroSect(band, Component_Particle)./ ExtCroSect(Band_Green, Component_Particle).* ...
            theta' / scale_factor;

        ssa_mixture = CompSSA(band, Component_Particle)*fraction_band';
        ssa_k = CompSSA(band, Component_Particle); %1x8
       
        tau = double(tau*scale_factor);

        tau_cam_k_ss = permute(reshape(smart.ss(:, ceil(x/RegScale),ceil(y/RegScale), Component_Particle, band, :),Model_OpticalDepthLen,Component_Num,Cam_Dim),[1,3,2]);
        tau_cam_k_ms = permute(reshape(smart.ms(:, ceil(x/RegScale),ceil(y/RegScale), Component_Particle, band, :),Model_OpticalDepthLen,Component_Num,Cam_Dim),[1,3,2]);
        rayleigh_k = reshape(smart.ms(1, ceil(x/RegScale),ceil(y/RegScale), Component_Particle, band, :),Component_Num,Cam_Dim)';

        ss = interpol_3d_mex(tau, tau_cam_k_ss); %9x8
        ms = interpol_3d_mex(tau, tau_cam_k_ms);

        atm_path(:, band) = (rayleigh_k + ss) * fraction_band' +  (ms - rayleigh_k) * (ssa_mixture ./ ssa_k .* exp(- tau * ...
            abs(ssa_mixture - ssa_k)) .* fraction_band)'; %9x1
    end

end