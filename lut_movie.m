function lut_movie(Date,Path,Orbit,Block,r,show_type,const)

        smart = load_cache(Date,Path,Orbit,Block,r,'smart');
        
        Orbit = num2str(Orbit,'%06d');
        Path = num2str(Path,'%03d');
        Block = num2str(Block);

        rho_ss = smart.ss;
        rho_ms = smart.ms;
        rho = rho_ss + rho_ms;
        cx = [0,0.7];

        %tmp = nanmean(nanmean(rho,3),2);
        tmp = rho(:,4,16,const.Component_Particle,2,:);
        rho_band = reshape(tmp,const.Model_OpticalDepthLen,8,const.Cam_Dim);
        tmp = rho(:,4,16,const.Component_Particle,:,1);
        rho_cam = reshape(tmp,const.Model_OpticalDepthLen,8,const.Band_Dim);
        
        cnt = 1;
        
        dir_video = 'video/';
        
        if ~exist(dir_video,'dir')
            mkdir(dir_video)
            fprintf('%s is created!\n',dir_video)
        else
            fprintf('directory %s exists, continue to save files!\n',dir_video)
        end
        
        file_video = strcat(Date,'_P',Path,'_O',Orbit,'_B',Block,'_',show_type,'.avi');
       
        writerObj = VideoWriter(strcat(dir_video,'/',file_video));
        writerObj.FrameRate = 3;
        open(writerObj)
        
        if strcmp(show_type,'cam_component')
            
            hFig = figure;
            set(hFig, 'Position', [300 300 1200 800])
        
            for aod_ind = 1:length(const.Model_OpticalDepthGrid)
                surf(squeeze(rho_band(aod_ind,:,:))),view([-150,30]),colorbar,caxis(cx)
                set(gca,'XTickLabel',const.Cam_Degree)
                set(gca,'YTickLabel',{'1','2','3','6','8','14','19','21'})
                zlim(cx)
                xlabel('Camera Off-Nadir Angle'),ylabel('Component Index')
                set(gca,'FontSize',18)
                title(strcat('\tau_{\lambda0}=',num2str(const.Model_OpticalDepthGrid(aod_ind))),'color','r','Fontsize',24)
                M = getframe(gcf);
                writeVideo(writerObj,M)
                cnt = cnt + 1;
            end
            
            
        elseif strcmp(show_type,'cam_aod')
            
            hFig = figure;
            set(hFig, 'Position', [300 300 1800 600])
            
            for i = 1:length(const.Component_Particle)
                surf(const.Model_OpticalDepthGrid,1:9,squeeze(rho_band(:,i,:))');view([-40,30]),colorbar,caxis(cx)
                zlim(cx)
                set(gca,'YTickLabel',const.Cam_Degree)
                ylabel('Camera Off-Nadir Angle'),xlabel('Aerosol Optical Depth')
                set(gca,'FontSize',18)
                title(strcat('Component:',num2str(const.Component_Particle(i))),'color','r','Fontsize',24)
                M = getframe(gcf);
                writeVideo(writerObj,M)
                cnt = cnt + 1;
            end
         
        elseif strcmp(show_type, 'band_component')
            
            hFig = figure;
            set(hFig, 'Position', [300 300 800 1000])
            
            for aod_ind = 1:length(const.Model_OpticalDepthGrid)
                surf([440,558,672,867],1:8,squeeze(rho_cam(aod_ind,:,:))),view([-150,30]),colorbar,caxis(cx)
                %set(gca,'XTickLabel',const.Band_Name)
                set(gca,'XTickLabel',{'Blue(440nm)','Green(558nm)','Red(672nm)','NIR(867nm)'})
                set(gca,'YTickLabel',{'1','2','3','6','8','14','19','21'})
                xlim([440,867])
                zlim(cx)
                xlabel('Band Name'),ylabel('Component Index')
                set(gca,'FontSize',18)
                title(strcat('\tau_{\lambda0}=',num2str(const.Model_OpticalDepthGrid(aod_ind))),'color','r','Fontsize',24)
                M = getframe(gcf);
                writeVideo(writerObj,M)
                cnt = cnt + 1;
            end
            
            
        else
            
            error('need to specify show type!\n')
                     
        end
                        
        close(writerObj)
        
end