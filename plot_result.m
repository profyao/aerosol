function plot_result(plot_name,r,const,varargin)
    
    cols = [0,0,1; 0,1,0; 1,0,0; 1,165/255,0];
    
    if ~exist('plots','dir')
        mkdir('plots')
        fprintf('directory plots is created!\n')
    end
    
    if strcmp(plot_name,'scatter_4band')
        Location = varargin{1};
        Method = varargin{2};
        
        if nargin~=5
            error('incorrect args for scatter_4band plot!\n')
        else
            figure
            for band = 1:const.Band_Dim
                aod_max = 0;
                subplot(2,2,band)
                [date_key,key,theta,aod_model,aod_aeronet] = load_aod_batch(Location,r,const,Method);
                aod_model_band = aod_model(:,band);
                aod_aeronet_band = aod_aeronet(:,band);
                aod_max = max([aod_max;aod_model_band;aod_aeronet_band]);
                h = scatter(aod_aeronet_band,aod_model_band,'MarkerEdgeColor',cols(band,:),'MarkerFaceColor',cols(band,:));
                fprintf('load band: %s\n', const.Band_Name{band})
                fprintf('correlation with aeronet is %f\n',corr(aod_model_band,aod_aeronet_band))
                fprintf('rms is %f\n',norm(aod_model_band-aod_aeronet_band)/sqrt(length(aod_model_band)))
                fprintf('avg bias is %f\n',mean(aod_model_band-aod_aeronet_band))
                aod_max_lim = 1.05*aod_max;
                xlim([0,aod_max_lim]),ylim([0,aod_max_lim])
                
                set(h,'SizeData',25)
                line('XData', [0 aod_max_lim], 'YData', [0 aod_max_lim], 'LineStyle', '-','LineWidth', 1, 'Color','k')

                xlabel('AERONET Measurement')
                ylabel('AOD Retrieval')
                title(strcat('Band:',const.Band_Name{band}))
                set(gca,'FontSize',18)
            end
        end
               
        %export_fig(strcat('plots/',plot_name,'_R',r,'_',strjoin(varargin,'_')),'-png','-transparent','-r240')
        
    elseif strcmp(plot_name,'scatter')
        
        Location = varargin{1};
        aod_max = 0;
        
        figure
        for p = 2:nargin - 3
            
            Method = varargin{p};          
            [date_key,key,xid,yid,lon,lat,theta,aod_model,aod_aeronet] = load_aod_batch(Location,r,const,Method);
    
            
            if ~strcmp(Method,'MISR')
                aod_model_band = aod_model(:,const.Band_Green);
            else
                aod_model_band = aod_model;
            end
            
            aod_aeronet_band = aod_aeronet(:,const.Band_Green);
            aod_max = max([aod_max;aod_model_band;aod_aeronet_band]);
            h = scatter(aod_aeronet_band,aod_model_band,'MarkerEdgeColor',const.cols(p-1,:),'MarkerFaceColor',const.cols(p-1,:));
            hold on
            fprintf('load band: %s\n', const.Band_Name{const.Band_Green})
            fprintf('correlation with aeronet is %f\n',corr(aod_model_band,aod_aeronet_band))
            fprintf('rms is %f\n',norm(aod_model_band-aod_aeronet_band)/sqrt(length(aod_model_band)))
            fprintf('avg bias is %f\n',mean(aod_model_band-aod_aeronet_band))
        end
        
        set(h,'SizeData',25)
        aod_max_lim = 1.05*aod_max;
        xlim([0,aod_max_lim]),ylim([0,aod_max_lim])
        line('XData', [0 aod_max_lim], 'YData', [0 aod_max_lim], 'LineStyle', '-','LineWidth', 1, 'Color','k')
        xlabel('AERONET Measurement')
        ylabel('AOD Retrieval')
        set(gca,'FontSize',18)
        legend(varargin{2:end},'Location','northwest')

    elseif strcmp(plot_name,'overlay')
        Date = varargin{1};
        Path = varargin{2};
        Orbit = varargin{3};
        Block = varargin{4};
        Location = varargin{5};
        Method = varargin{6};
        cmap = varargin{7};
        
        [aod_model, ~,~, ~, lon1,lat1] = load_aod(Date,Path,Orbit,Block,r,const,Method);
        [aod_aeronet, ~, ~, lon2,lat2] = load_aeronet(Date,Path,Block,r,Location,const);
        
        if ~strcmp(Method,'MISR')
            aod_model = aod_model(:,2);
        end
        aod_aeronet = aod_aeronet(:,2);
        
        %aod = [aod_model;aod_aeronet];
        %aod_min = min(aod);aod_max = max(aod);
        aod_min = 0.02;aod_max=0.2;
        figure
        scatter(lon1,lat1,100,aod_model,'s','filled'), set(gca,'CLim',[aod_min aod_max]), hold on
        scatter(lon2,lat2,100,aod_aeronet,'o','filled','MarkerEdgeColor','k'), set(gca,'CLim',[aod_min aod_max])
        colormap(cmap),colorbar
        %uistack(h,'bottom');
        
        plot_google_map('MapType','hybrid','ShowLabels',0,'Alpha',0.8)
        set(gca,'FontSize',18)
        file_model = strcat('plots/',plot_name,'_R',r,strjoin({Date,Location,Method},'_'));
        %export_fig(file_model,'-png','-transparent','-r240')
    
    elseif strcmp(plot_name,'resid')
        
        Date = varargin{1};
        Path = varargin{2};
        Orbit = varargin{3};
        Block = varargin{4};
        Method = varargin{5};
        iter = varargin{6};
        
        sample = load_cache(Date,Path,Orbit,Block,r,'sample',Method);
        resid = squeeze(sample.resid(:,:,iter))';
        resid(isinf(resid)) = NaN;
       
        C = nancov(resid);
        [R,~] = corrcov(C);
        imagesc(R),colorbar;
        set(gcf,'Position',[100, 100, 800, 600]);
        set(gca,'FontSize',18)
        %export_fig(strcat('plots/',plot_name,'_R',rstrjoin({Date,Method,num2str(iter)})),'-png','-transparent','-r240')
        
    elseif strcmp(plot_name,'stability')
       
        Location = varargin{1};
        test_delta = varargin{2};
        
        figure
        aod_max = 0;
        
        for p = 3:length(varargin)
            Method = varargin{p};

            [aoda,aodb] = load_aod_rep(Location,r,const,Method,test_delta);
            
            if test_delta == 0
                avga  = mean(aoda,2);
                %stda = std(aoda,0,2);
                bd = quantile(aoda,[0.1,0.9],2);
                errorbar(aodb,avga,bd(:,1)-avga,bd(:,2)-avga,'MarkerEdgeColor',const.cols(p-2,:),'MarkerFaceColor',const.cols(p-2,:),'Color',const.cols(p-2,:),'LineStyle','none','Marker','x')
                aod_max = max([aod_max;bd(:,2);aodb]); 
            else
                plot(repmat(aodb,1,6),aoda,'x')
                legend({'\Delta=0.001','\Delta=0.01','\Delta=0.05','\Delta=0.1','\Delta=1','\Delta=10'},'Location','southeast')
                aod_max = max([aod_max;aoda(:);aodb]);
                delta_all = [0.001,0.01,0.05,0.1,1,10];
                
                for rep = 1:6
                    fprintf('delta is %f\n',delta_all(rep))
                    fprintf('correlation with aeronet is %f\n',corr(aodb,aoda(:,rep)))
                    fprintf('rms is %f\n',norm(aodb-aoda(:,rep))/sqrt(length(aodb)))
                    fprintf('avg bias is %f\n',mean(aodb-aoda(:,rep)))
                end
                
            end
            
            hold on
        end
        aod_max_lim = 1.05*aod_max;
        xlim([0,aod_max_lim]),ylim([0,aod_max_lim])
        line('XData', [0 aod_max_lim], 'YData', [0 aod_max_lim], 'LineStyle', '-','LineWidth', 1, 'Color','k')
        if test_delta == 0
            legend({'Random Local Search','MCMC'},'Location','northwest')
        end
        xlabel('AERONET Measurement','FontSize',18)
        ylabel('AOD Retrieval','FontSize',18)
        set(gca,'FontSize',18)
        if test_delta == 0
            file_model = strcat('plots/',plot_name,'_boot');
        else
            file_model = strcat('plots/',plot_name,'_delta');
        end
        %export_fig(file_model,'-png','-transparent','-r240')
        
        
    elseif strcmp(plot_name,'post_dist')
        
        Location = varargin{1};
        idx = varargin{2}; % 1 to 88
        [aoda_mcmc,aodb] = load_aod_rep(Location,const,'MCMC',0);
        [aoda_cdr,~] = load_aod_rep(Location,const,'CD-random',0);
        
        file_video = strcat('video/post_dist.avi');
        writerObj = VideoWriter(file_video);
        writerObj.FrameRate = 3;
        open(writerObj)
        
        figure
        for id = 1:88
            xrange = [min([aoda_mcmc(id,:),aoda_cdr(id,:),aodb(id)]),max([aoda_mcmc(id,:),aoda_cdr(id,:),aodb(id)])];
            subplot(2,1,1)
            hist(aoda_mcmc(id,:))
            xlim(xrange)
            vline(aodb(id),'r')
            title('MCMC Retrieval','FontSize',18)
            ylabel('Frequency','FontSize',18)
            xlabel(strcat('Pixel Number:',num2str(id)),'FontSize',18)
            subplot(2,1,2)
            hist(aoda_cdr(id,:))
            xlim(xrange)
            vline(aodb(id),'r')
            title('Random Local Search Retrieval','FontSize',18)
            ylabel('Frequency','FontSize',18)
            xlabel(strcat('Pixel Number:',num2str(id)),'FontSize',18)

            
            M = getframe(gcf);
            writeVideo(writerObj,M)
        end
        close(writerObj)
        
    elseif strcmp(plot_name,'theta')
        
        Date = varargin{1};
        Path = varargin{2};
        Orbit = varargin{3};
        Block = varargin{4};
        Method = varargin{5};
        cmap = varargin{6};
        
        [~,theta,xid,yid,~,~] = load_aod(Date,Path,Orbit,Block,r,const,Method);
        figure

        if ~strcmp(Method,'MISR')

            for j = 1:const.Component_Num
                h = subplot(4,2,j);
                p = get(h, 'pos');
                p(1) = p(1) - 0.05;
                p(3) = p(3) + 0.05;
                set(h, 'pos', p);

                plot_1d(r, theta(j,:), xid, yid, cmap, const,[0,0.5])

                title(strcat('Component Num:',num2str(const.Component_Particle(j))))
                set(gca,'FontSize',18)
            end 

        else

            for j = 1:5
                if j<=3
                    h = subplot(3,2,j);
                else
                    h = subplot(3,2,j+1);
                end
                p = get(h, 'pos');
                p(1) = p(1) - 0.05;
                p(3) = p(3) + 0.05;
                set(h, 'pos', p);

                plot_1d(r,theta(j,:),xid,yid,cmap,const,[0,1])

                switch j
                    case 1
                        title('Small Particle (<0.35\mu)')
                    case 2
                        title('Medium Particle (0.35\mu~0.7\mu)')
                    case 3
                        title('Large Particle (>0.7\mu)')
                    case 4
                        title('Spherical')
                    case 5
                        title('Non-Spherical')
                end

                set(gca,'FontSize',18)
            end

        end
    
        
    elseif strcmp(plot_name, 'trace')
        
        Location = varargin{1};
        Method = varargin{2};
        [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
        id = find(strcmp(TXT(2:end,7),Location));

        Dates = TXT(id+1,2);
        Paths = NUM(id+1,3);
        Orbits = NUM(id+1,4);
        Blocks = NUM(id+1,5);

        N = length(Dates);
        figure
        
        for i = 1:N
           
            Date = Dates{i};
            Path = Paths(i);
            Orbit = Orbits(i);
            Block = Blocks(i);
            [reg,sample] = load_cache(Date,Path,Orbit,Block,r,'reg','sample','MCMC');
            [x1,y1] = find(reg.reg_is_used);
            [~, x2, y2, lon_a, lat_a] = load_aeronet(Date,Path,Block,r,Location,const);
            [~,~,I1,I2] = match_aeronet(x1,y1,x2,y2);
            fprintf('MCMC:%d,aeronet:%d,%d points are matched!\n',length(x1),length(x2),length(I1))
            aod = sample.tau(I1,:);
            for j = 1:length(I1)
                if strcmp(Method,'mean')
                    plot(cumsum(aod(j,:))./[1:1001])
                elseif strcmp(Method,'raw')
                    tmp = aod(j,:);
                    dif = arrayfun(@(x) tmp(x) - tmp(x-1), 2:length(tmp));
                    disp(sum(dif~=0)/(length(tmp)-1))
                    plot(tmp)
                elseif strcmp(Method,'acf')
                    plot(acf(aod(j,:)',1000))
                end
                hold on
            end
        end
        
    end

end