function polar_movie(aod_ind,sfc_pres_ind,band_ind,scatter_type,const)
        
        dir_video = 'video/';
        
        if ~exist(dir_video,'dir')
            mkdir(dir_video);
            fprintf('%s is created!\n',dir_video)
        else
            fprintf('directory %s exists, continue to save files!\n',dir_video)
        end
        
        file_video = strcat('polar_',scatter_type,'.avi');
        
        writerObj = VideoWriter(strcat(dir_video,'/',file_video));
        open(writerObj);
        figure
        for component_ind = 1:length(const.Component_Particle)
            for mu0_ind = 1:length(const.Model_mu0Grid)
                clf;
                theta0 = polar_plot(aod_ind,sfc_pres_ind,band_ind,component_ind,mu0_ind,scatter_type,const);
                M = getframe(gcf);
                writeVideo(writerObj,M);
            end
        end
        
        close(writerObj)

end