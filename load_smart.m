function smart = load_smart(Date,Path,Orbit,Block,const,add_limit)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');
    
    file_smart_ss = fullfile('products/MIANSMT/',const.MIANSMT_SS_filename);
    file_smart_ms = fullfile('products/MIANSMT/',const.MIANSMT_MS_filename);

    AlgTypeFlag = hdfread(file_aerosol, 'RegParamsAer', 'Fields', 'AlgTypeFlag', ...
        'Index',{[Block  1  1],[1  1  1],[1  const.XDim_r17600  const.YDim_r17600]}); %8x32, 17.6km
    SolZenAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'SolZenAng', ...
    'Index',{[Block  1  1],[1  1  1],[1  const.XDim_r17600  const.YDim_r17600]}); % 8x32, 17.6km
    ViewZenAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'ViewZenAng', ...
        'Index',{[Block  1  1  1],[1  1  1  1],[1  const.XDim_r17600  ...
        const.YDim_r17600  const.Cam_Dim]}); % 8x32x9, 17.6km
    ScatterAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'ScatterAng', ...
        'Index',{[Block  1  1  1],[1  1  1  1],[1  const.XDim_r17600 ...
        const.YDim_r17600  const.Cam_Dim]}); % 8x32x9, 17.6km
    SfcPres = hdfread(file_aerosol, 'RegParamsEnvironmental', 'Fields', 'SfcPres', ...
        'Index',{[Block  1  1],[1  1  1],[1   const.XDim_r17600  const.YDim_r17600]}); %8x32
    
    if add_limit == true
        
        file_smart_tdiff = fullfile('products/MIANSMT/',const.MIANSMT_TDIFF_filename);
        file_smart_ediff = fullfile('products/MIANSMT/',const.MIANSMT_EDIFF_filename);
        
        EdiffSS = hdfread(file_smart_ediff, '/EdiffSingleScatter', 'Index', {[1  1  1  1  1],[1  1  1  1  1],[4  21  13   2  81]});
        EdiffMS = hdfread(file_smart_ediff, '/EdiffMultipleScatter', 'Index', {[1  1  1  1  1],[1  1  1  1  1],[4  21  13   2  81]});
        EdiffSS = exp(double(EdiffSS/1000));
        EdiffMS = exp(double(EdiffMS/1000));
        TdiffSS = hdfread(file_smart_tdiff,'/TSingleScatter', 'Index', {[1  1  1  1  1],[1  1  1  1  1],[4  21  13   2  29]});
        TdiffMS = hdfread(file_smart_tdiff, '/TMultipleScatter', 'Index', {[1  1  1  1  1],[1  1  1  1  1],[4  21  13   2  29]});            
        TdiffSS = exp(double(TdiffSS/1000));
        TdiffMS = exp(double(TdiffMS/1000));
        
        smart.ss_tdiff = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim, const.Cam_Dim);
        smart.ss_ediff = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim);
        smart.ms_tdiff = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim, const.Cam_Dim);
        smart.ms_ediff = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim);

    end
    
    % interp view and scatter angles
    smart.ss = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim, const.Cam_Dim);
    smart.ms = NaN*ones(const.Model_OpticalDepthLen, const.XDim_r17600,const.YDim_r17600, const.Model_ComponentDim, const.Band_Dim, const.Cam_Dim);
    
    smart.mu = NaN*ones(const.XDim_r17600, const.YDim_r17600, const.Cam_Dim); % 8x32x9
    smart.scatter_angle = NaN*ones(const.XDim_r17600, const.YDim_r17600, const.Cam_Dim); % 8x32x9
    smart.mu0 = NaN*ones(const.XDim_r17600,const.YDim_r17600); %8x32
    
    fprintf('load smart data!\n')
    
    for ii = 1:const.XDim_r17600
        for jj = 1:const.YDim_r17600
            
            % view angle
            mu = NaN * ones(const.Cam_Dim,1);
            mu_ind = NaN * ones(const.Cam_Dim,1);
            % scatter angle
            scatter_angle = NaN * ones(const.Cam_Dim,1);
            scatter_angle_ind = NaN * ones(const.Cam_Dim,1);
            % pressure: trim outliers
            SfcPres(SfcPres > const.Model_Pressure(2)) = const.Model_Pressure(2);
            SfcPres(SfcPres < const.Model_Pressure(1)) = const.Model_Pressure(1);

            if AlgTypeFlag(ii,jj) == 3 %heterogeneous surface retrieval

                %disp([ii,jj]);

                mu0 = double(cosd(SolZenAng(ii,jj)));
                mu0_ind = find(const.Model_mu0Grid >= mu0, 1, 'first');
                surface_pressure = double(SfcPres(ii,jj));  

                % Load top-of-atmosphere equivalent reflectance
                data1 = hdfread(file_smart_ss,  ['/EquivalentReflectanceMu0_',num2str(mu0_ind-1)], 'Index', ...
                    {[1  1  1  1  1  1],[1  1  1  1  1  1],[4  21  13  2  29  96]});
                data2 = hdfread(file_smart_ss,  ['/EquivalentReflectanceMu0_',num2str(mu0_ind)], 'Index', ...
                    {[1  1  1  1  1  1],[1  1  1  1  1  1],[4  21  13  2  29  96]});
                SS = cat(ndims(data1)+1, data1, data2);

                data1 = hdfread(file_smart_ms,  ['/EquivalentReflectanceMu0_',num2str(mu0_ind-1)], 'Index', ...
                    {[1  1  1  1  1  1],[1  1  1  1  1  1],[4  21  13  2  29  96]});
                data2 = hdfread(file_smart_ms,  ['/EquivalentReflectanceMu0_',num2str(mu0_ind)], 'Index', ...
                    {[1  1  1  1  1  1],[1  1  1  1  1  1],[4  21  13  2  29  96]});
                MS = cat(ndims(data1)+1, data1, data2);
                
                % build up data set for different cams
                ss_cam = NaN * ones(const.Cam_Dim, const.Band_Dim, const.Model_ComponentDim, const.Model_OpticalDepthLen,...
                    length(const.Model_Pressure), 2, 4, 2);
                ms_cam = NaN * ones(const.Cam_Dim, const.Band_Dim, const.Model_ComponentDim, const.Model_OpticalDepthLen,...
                    length(const.Model_Pressure), 2, 4, 2);

                for cam = 1:const.Cam_Dim

                    mu(cam) = double(cosd(ViewZenAng(ii, jj, cam)));
                    mu_ind(cam) = find(const.Model_muGrid>=mu(cam), 1, 'first');

                    scatter_angle(cam) = double(ScatterAng(ii,jj,cam));
                    scatter_angle_ind(cam) = find(const.Model_ScatterAngleGrid>=scatter_angle(cam), 1, 'first');

                    ss_cam(cam,:,:,:,:,:,:,:) = SS(:,:,:,:,[mu_ind(cam)-1, mu_ind(cam)],[1, scatter_angle_ind(cam)-1, scatter_angle_ind(cam), end], :);
                    ms_cam(cam,:,:,:,:,:,:,:) = MS(:,:,:,:,[mu_ind(cam)-1, mu_ind(cam)],[1, scatter_angle_ind(cam)-1, scatter_angle_ind(cam), end], :);

                    I = find(reshape(ss_cam(cam,1,1,1,1,:,:,:),2,4,2)==32767|reshape(ms_cam(cam,1,1,1,1,:,:,:),2,4,2)==32767);

                    if ~isempty(I)
                        [J,K,L] = ind2sub([2,4,2],I);
                        low_allowable = 180 - abs(ViewZenAng(ii,jj,cam)+SolZenAng(ii,jj));
                        high_allowable = 180 - abs(ViewZenAng(ii,jj,cam)-SolZenAng(ii,jj));

                        if scatter_angle(cam) > (low_allowable + high_allowable)/2
                            for kk = 1:length(I)
                                ss_cam(cam,:,:,:,:,J(kk),K(kk),L(kk)) = ss_cam(cam,:,:,:,:,J(kk),end,L(kk));
                                ms_cam(cam,:,:,:,:,J(kk),K(kk),L(kk)) = ms_cam(cam,:,:,:,:,J(kk),end,L(kk));
                            end
                        else
                            for kk = 1:length(I)
                                ss_cam(cam,:,:,:,:,J(kk),K(kk),L(kk)) = ss_cam(cam,:,:,:,:,J(kk),1,L(kk));
                                ms_cam(cam,:,:,:,:,J(kk),K(kk),L(kk)) = ms_cam(cam,:,:,:,:,J(kk),1,L(kk));
                            end
                        end
                    end
                end

                ss_cam = exp(double(ss_cam)/1000);
                ms_cam = exp(double(ms_cam)/1000);

                for cam = 1:const.Cam_Dim
                    [x1,x2,x3,x4,x5] = ndgrid(const.Model_OpticalDepthGrid, const.Model_Pressure,...
                                const.Model_muGrid([mu_ind(cam)-1, mu_ind(cam)]), ...
                                const.Model_ScatterAngleGrid([scatter_angle_ind(cam)-1, scatter_angle_ind(cam)]), ...
                                const.Model_mu0Grid([mu0_ind-1, mu0_ind]));

                    for kk = 1:const.Model_ComponentDim
                        for band  = 1:const.Band_Dim %opt_depth x pressure x mu x 
                            smart.ss(:, ii,jj,kk, band, cam) = interpn(x1,x2,x3,x4,x5,reshape(ss_cam(cam, band, kk,:,:,:,[2,3],:),const.Model_OpticalDepthLen,2,2,2,2),...
                                const.Model_OpticalDepthGrid, surface_pressure, mu(cam), scatter_angle(cam), mu0);
                            smart.ms(:, ii,jj,kk, band, cam) = interpn(x1,x2,x3,x4,x5,reshape(ms_cam(cam, band, kk,:,:,:,[2,3],:),const.Model_OpticalDepthLen,2,2,2,2),...
                                const.Model_OpticalDepthGrid, surface_pressure, mu(cam), scatter_angle(cam), mu0);   
                        end
                    end
                end
                
                if add_limit == true
                    % Load diffuse irradiance at the bottom of the atmosphere
                    SS = EdiffSS(:,:,:,:,[mu0_ind-1,mu0_ind]);
                    MS = EdiffMS(:,:,:,:,[mu0_ind-1,mu0_ind]);
                    [x1,x2,x3] = ndgrid(const.Model_OpticalDepthGrid, const.Model_Pressure, const.Model_mu0Grid([mu0_ind-1,mu0_ind]));
                    for kk = 1:const.Model_ComponentDim
                        for band  = 1:const.Band_Dim %opt_depth x pressure x mu x 
                            smart.ss_ediff(:, ii,jj, kk, band) = interpn(x1,x2,x3,squeeze(SS(band, kk,:,:,:)),...
                                const.Model_OpticalDepthGrid, surface_pressure, mu0);
                            smart.ms_ediff(:, ii,jj, kk, band) = interpn(x1,x2,x3,squeeze(MS(band, kk,:,:,:)),...
                                const.Model_OpticalDepthGrid, surface_pressure, mu0);
                        end
                    end

                    % build up data set for different cams

                    ss_cam = NaN * ones(const.Cam_Dim, const.Band_Dim, const.Model_ComponentDim, const.Model_OpticalDepthLen,...
                    length(const.Model_Pressure), 2);
                    ms_cam = NaN * ones(const.Cam_Dim, const.Band_Dim, const.Model_ComponentDim, const.Model_OpticalDepthLen,...
                    length(const.Model_Pressure), 2);

                    for cam = 1:const.Cam_Dim
                        ss_cam(cam,:,:,:,:,:) = TdiffSS(:,:,:,:,[mu_ind(cam)-1,mu_ind(cam)]);
                        ms_cam(cam,:,:,:,:,:) = TdiffMS(:,:,:,:,[mu_ind(cam)-1,mu_ind(cam)]);      
                    end

                    for cam = 1:const.Cam_Dim    
                        [x1,x2,x3] = ndgrid(const.Model_OpticalDepthGrid, const.Model_Pressure, const.Model_muGrid([mu_ind(cam)-1,mu_ind(cam)]));
                        for kk = 1:const.Model_ComponentDim
                            for band = 1:const.Band_Dim
                                smart.ss_tdiff(:, ii,jj,kk, band, cam)=interpn(x1,x2,x3, squeeze(ss_cam(cam,band,kk,:,:,:)),...
                                    const.Model_OpticalDepthGrid, surface_pressure, mu(cam));
                                smart.ms_tdiff(:, ii,jj,kk, band, cam)=interpn(x1,x2,x3, squeeze(ms_cam(cam,band,kk,:,:,:)),...
                                    const.Model_OpticalDepthGrid, surface_pressure, mu(cam));

                            end
                        end

                    end
                end

                smart.mu0(ii,jj) = mu0;
                smart.mu(ii,jj,:) = mu;
                smart.scatter_angle(ii,jj,:) = scatter_angle;
            else
                %fprintf('%d,%d: not heterogeneous region!\n',ii,jj)
            end    
        end
    end 
end