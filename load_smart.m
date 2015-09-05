function smart = load_smart(Date,Path,Orbit,Block,const)
    
    MIANSMT_MS_filename = const.MIANSMT_MS_filename;
    MIANSMT_SS_filename = const.MIANSMT_SS_filename;
    header_MIL2ASAE_filename = const.header_MIL2ASAE_filename;
    Cam_Dim = const.Cam_Dim;
    Band_Dim = const.Band_Dim;
    Model_ComponentDim = const.Model_ComponentDim;
    Model_OpticalDepthLen = const.Model_OpticalDepthLen;
    Model_OpticalDepthGrid = const.Model_OpticalDepthGrid;
    Model_mu0Grid = const.Model_mu0Grid;
    Model_muGrid = const.Model_muGrid;
    Model_ScatterAngleGrid = const.Model_ScatterAngleGrid;
    Model_Pressure = const.Model_Pressure;
    XDim_r17600 = const.XDim_r17600;
    YDim_r17600 = const.YDim_r17600;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');
    
    file_smart_ss = fullfile('products/MIANSMT/',MIANSMT_SS_filename);
    file_smart_ms = fullfile('products/MIANSMT/',MIANSMT_MS_filename);

    AlgTypeFlag = hdfread(file_aerosol, 'RegParamsAer', 'Fields', 'AlgTypeFlag', ...
        'Index',{[Block  1  1],[1  1  1],[1  XDim_r17600  YDim_r17600]}); %8x32, 17.6km
    SolZenAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'SolZenAng', ...
    'Index',{[Block  1  1],[1  1  1],[1  XDim_r17600  YDim_r17600]}); % 8x32, 17.6km
    ViewZenAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'ViewZenAng', ...
        'Index',{[Block  1  1  1],[1  1  1  1],[1  XDim_r17600  ...
        YDim_r17600  Cam_Dim]}); % 8x32x9, 17.6km
    ScatterAng = hdfread(file_aerosol, 'RegParamsGeometry', 'Fields', 'ScatterAng', ...
        'Index',{[Block  1  1  1],[1  1  1  1],[1  XDim_r17600 ...
        YDim_r17600  Cam_Dim]}); % 8x32x9, 17.6km
    SfcPres = hdfread(file_aerosol, 'RegParamsEnvironmental', 'Fields', 'SfcPres', ...
        'Index',{[Block  1  1],[1  1  1],[1   XDim_r17600  YDim_r17600]}); %8x32
    
    % build up data set for different cams
    ss_cam = NaN * ones(Cam_Dim, Band_Dim, Model_ComponentDim, Model_OpticalDepthLen,...
        length(Model_Pressure), 2, 4, 2);
    ms_cam = NaN * ones(Cam_Dim, Band_Dim, Model_ComponentDim, Model_OpticalDepthLen,...
        length(Model_Pressure), 2, 4, 2);
    
    % interp view and scatter angles
    smart.ss = NaN*ones(Model_OpticalDepthLen, XDim_r17600,YDim_r17600, Model_ComponentDim, Band_Dim, Cam_Dim);
    smart.ms = NaN*ones(Model_OpticalDepthLen, XDim_r17600,YDim_r17600, Model_ComponentDim, Band_Dim, Cam_Dim);
    
    smart.mu = NaN*ones(XDim_r17600, YDim_r17600, Cam_Dim); % 8x32x9
    smart.scatter_angle = NaN*ones(XDim_r17600, YDim_r17600, Cam_Dim); % 8x32x9
    smart.mu0 = NaN*ones(XDim_r17600,YDim_r17600); %8x32
    
    fprintf('load smart data!\n')
    
    for ii = 1:XDim_r17600
        for jj = 1:YDim_r17600
            
            % view angle
            mu = NaN * ones(Cam_Dim,1);
            mu_ind = NaN * ones(Cam_Dim,1);
            % scatter angle
            scatter_angle = NaN * ones(Cam_Dim,1);
            scatter_angle_ind = NaN * ones(Cam_Dim,1);
            % pressure: trim outliers
            SfcPres(SfcPres > Model_Pressure(2)) = Model_Pressure(2);
            SfcPres(SfcPres < Model_Pressure(1)) = Model_Pressure(1);

            if AlgTypeFlag(ii,jj) == 3 %heterogeneous surface retrieval

                %disp([ii,jj]);

                mu0 = double(cosd(SolZenAng(ii,jj)));
                mu0_ind = find(Model_mu0Grid >= mu0, 1, 'first');
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

                for cam = 1:Cam_Dim

                    mu(cam) = double(cosd(ViewZenAng(ii, jj, cam)));
                    mu_ind(cam) = find(Model_muGrid>=mu(cam), 1, 'first');

                    scatter_angle(cam) = double(ScatterAng(ii,jj,cam));
                    scatter_angle_ind(cam) = find(Model_ScatterAngleGrid>=scatter_angle(cam), 1, 'first');

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

                for cam = 1:Cam_Dim
                    [x1,x2,x3,x4,x5] = ndgrid(Model_OpticalDepthGrid, Model_Pressure,...
                                Model_muGrid([mu_ind(cam)-1, mu_ind(cam)]), ...
                                Model_ScatterAngleGrid([scatter_angle_ind(cam)-1, scatter_angle_ind(cam)]), ...
                                Model_mu0Grid([mu0_ind-1, mu0_ind]));

                    for kk = 1:Model_ComponentDim
                        for band  = 1:Band_Dim %opt_depth x pressure x mu x 
                            smart.ss(:, ii,jj,kk, band, cam) = interpn(x1,x2,x3,x4,x5,reshape(ss_cam(cam, band, kk,:,:,:,[2,3],:),Model_OpticalDepthLen,2,2,2,2),...
                                Model_OpticalDepthGrid, surface_pressure, mu(cam), scatter_angle(cam), mu0);
                            smart.ms(:, ii,jj,kk, band, cam) = interpn(x1,x2,x3,x4,x5,reshape(ms_cam(cam, band, kk,:,:,:,[2,3],:),Model_OpticalDepthLen,2,2,2,2),...
                                Model_OpticalDepthGrid, surface_pressure, mu(cam), scatter_angle(cam), mu0);   
                        end
                    end
                end
                
                smart.mu0(ii,jj) = mu0;
                smart.mu(ii,jj,:) = mu;
                smart.scatter_angle(ii,jj,:) = scatter_angle;
            else
                fprintf('%d,%d: not heterogeneous region!\n',ii,jj)
            end    
        end
    end 
end