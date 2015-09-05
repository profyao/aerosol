function theta0 = polar_plot(aod_ind,srf_pres_ind,band_ind,component_ind,mu0_ind,scatter_type,const)
    
    [theta0,theta,~,phi,rho] = extract_BRF(aod_ind,srf_pres_ind,band_ind,component_ind,mu0_ind,scatter_type,const);
    aod_grid = const.Model_OpticalDepthGrid;
    
    % Convert to Cartesian
    %x = theta.*cosd(phi);
    %y = theta.*sind(phi);
    
    x = sind(theta).*cosd(phi);
    y = sind(theta).*sind(phi);
    z = cosd(theta);
    
    scatter3([x;x],[y;-y],[z;z],[],[rho;rho],'filled');colorbar;view([10,30]),caxis([0,0.2])
    hold on
    plot3(sin(pi/2:-pi/180:-pi/2),zeros(1,181),cos(pi/2:-pi/180:-pi/2),'--r');
    hold on
    plot3(zeros(1,181),sin(pi/2:-pi/180:-pi/2),cos(pi/2:-pi/180:-pi/2),'--r');
    hold on
    plot3([-1,1],[0,0],[0,0],'--r');
    hold on
    plot3([0,0],[-1,1],[0,0],'--r');
    hold on
    plot3([0,0],[0,0],[0,1],'r','LineWidth',4);
    %h = polar(phi*pi/180,theta);
    hold on
    plot3(cos(0:2*pi/360:2*pi),sin(0:2*pi/360:2*pi),zeros(361,1),'--r');
    hold on
    scatter3(sind(theta0),0,cosd(theta0),'+r');
    hold on
    mArrow3([sind(theta0),0,cosd(theta0)],[0,0,0], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 0.01);
    text(1.1*sind(theta0),0,1.1*cosd(theta0),num2str(theta0,3),'color','r','Fontsize',18);
    %contourf(x,y,z),colorbar;
    %scatter([x;x],[y;-y],[],[rho;rho],'.'),colorbar;
    %hold on
    %polar(0,theta0,'+r')
    title(strcat('Component Particle:',num2str(component_ind),'Band:',num2str(band_ind),'AOD:',num2str(aod_grid(aod_ind))),...
        'color','r','Fontsize',18);
    % Hide the POLAR function data and leave annotations
    %set(h,'Visible','off')
    % Turn off axes and set square aspect ratio
    %axis off
    %axis image
    %text(theta0,0,num2str(theta0,3),'color','r','Fontsize',18)

end