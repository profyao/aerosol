function lut_movie(Date,Path,Orbit,Block,r,scatter_type,show_type,const)

        smart = load_cache(Date,Path,Orbit,Block,r,'smart');
        
        Orbit = num2str(Orbit,'%06d');
        Path = num2str(Path,'%03d');
        Block = num2str(Block);

        if strcmp(scatter_type,'SS')
            rho = smart.ss;
            cx = [0,0.15];
        elseif strcmp(scatter_type,'MS')
            rho = smart.ms;
            cx = [0,0.7];
        else
            error('need to specify scattering type!\n')
        end
        
        %tmp = nanmean(nanmean(rho,3),2);
        tmp = rho(:,4,16,:,:,:);
        rho = reshape(tmp,const.Model_OpticalDepthLen,const.Model_ComponentDim,const.Band_Dim,const.Cam_Dim);
        
        cnt = 1;
        
        dir_video = 'video/';
        
        if ~exist(dir_video,'dir')
            mkdir(dir_video)
            fprintf('%s is created!\n',dir_video)
        else
            fprintf('directory %s exists, continue to save files!\n',dir_video)
        end
        
        file_video = strcat(Date,'_P',Path,'_O',Orbit,'_B',Block,'_',show_type,'_',scatter_type,'.avi');
       
        writerObj = VideoWriter(strcat(dir_video,'/',file_video));
        writerObj.FrameRate = 3;
        open(writerObj)
        
        figure
        
        if strcmp(show_type,'aod_component')
       
            for band = 1:const.Band_Dim
                for cam = 1:const.Cam_Dim
                    surf(rho(:,:,band,cam)'),view([-80,30]),colorbar,caxis(cx)
                    title(strcat('Aerosol Reflectance:','Band',num2str(band),'Cam',num2str(cam)),'color','r','Fontsize',18)
                    set(gca,'XTickLabel',strtrim(cellstr(num2str(const.Model_OpticalDepthGrid(2:3:13)'))'))
                    xlabel('Aerosol Optical Depth'),ylabel('Component Index')
                    M = getframe(gcf);
                    writeVideo(writerObj,M)
                    cnt = cnt + 1;
                end
            end
            
        elseif strcmp(show_type,'aod_cam')
            
            for band = 1:const.Band_Dim
                for component_ind = const.Component_Particle
                    surf(reshape(rho(:,component_ind,band,:),const.Model_OpticalDepthLen,const.Cam_Dim)');view([-80,30]),colorbar,caxis(cx)
                    zlim(cx)
                    title(strcat('Aerosol Reflectance:','Band',num2str(band),'Component',num2str(component_ind)),'color','r','Fontsize',18)
                    set(gca,'XTickLabel',strtrim(cellstr(num2str(const.Model_OpticalDepthGrid(2:3:13)'))'))
                    xlabel('Aerosol Optical Depth'),ylabel('Cam Index')
                    M = getframe(gcf);
                    writeVideo(writerObj,M)
                    cnt = cnt + 1;
                end
            end
            
        else
            
            error('need to specify show type!\n')
                     
        end
                        
        close(writerObj)
        
end