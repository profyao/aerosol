function aod = extract_aod_4band(tau,theta,num_reg_used,ExtCroSect,const)

    aod = zeros(num_reg_used,const.Band_Dim);
    
    for p = 1:num_reg_used
        for band = 1:const.Band_Dim
            scale_factor = ExtCroSect(band, const.Component_Particle) ./ ExtCroSect(const.Band_Green, const.Component_Particle) * theta(:,p);
            aod(p,band) = tau(p)*scale_factor;
        end
    end

end