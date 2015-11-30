function get_image(Date,Path,Orbit,Block,r,Band,const,simulation)
    
    if simulation == 1
        reg_sim = load_cache(Date,Path,Orbit,Block,r,'reg_sim','sim');
        equ_ref = reg_sim;
    else
        reg = load_cache(Date,Path,Orbit,Block,r,'reg');
        equ_ref = reg.mean_equ_ref;
    end
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_video = 'video/';
        
    if ~exist(dir_video,'dir')
        mkdir(dir_video)
        fprintf('%s is created!\n',dir_video)
    else
        fprintf('directory %s exists, continue to save files!\n',dir_video)
    end
    
    if simulation == 0
        file_video = strcat(Date,'_P',Path,'_O',Orbit,'_B',Block,'_image.avi');
    else
        file_video = strcat(Date,'_P',Path,'_O',Orbit,'_B',Block,'_image_sim.avi');
    end

    writerObj = VideoWriter(strcat(dir_video,'/',file_video));
    writerObj.FrameRate = 1;
    open(writerObj)
    figure('position',[300 300 1200 300])
    %cx = [min(equ_ref(:)),max(equ_ref(:))];
    cx = [0.06,0.18];
    
    for cam = 1:const.Cam_Dim
        plot_2d(equ_ref(:,:,Band,cam),jet,cx)
        title(sprintf('Band: %s, Cam: %d',const.Band_Name{Band},cam),'Fontsize',24)
        M = getframe(gcf);
        writeVideo(writerObj,M)
    end
    
    close(writerObj)

end