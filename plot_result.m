function plot_result(plot_name,const,varargin)
    
    cols = const.cols;
    
    if ~exist('plots','dir')
        mkdir('plots')
        fprintf('directory plots is created!\n')
    end
    
    if strcmp(plot_name,'scatter')
        Location = varargin{1};
        aod_max = 0;
        
        if nargin<4
            error('not enough args for scatter plot!\n')
        else
            figure
            for p = 2:length(varargin)
                [aod_model,aod_aeronet] = load_aod_batch(Location,const,varargin{p});
                aod_max = max([aod_max;aod_model;aod_aeronet]);
                h = scatter(aod_aeronet,aod_model,'MarkerEdgeColor',cols(p-1,:),'MarkerFaceColor',cols(p-1,:));
                fprintf('correlation with aeronet is %f\n',corr(aod_model,aod_aeronet))
                fprintf('rms is %f\n',norm(aod_model-aod_aeronet)/sqrt(length(aod_model)))
                fprintf('avg bias is %f\n',mean(aod_model-aod_aeronet))
                hold on
            end
            aod_max_lim = 1.05*aod_max;
            xlim([0,aod_max_lim]),ylim([0,aod_max_lim])
            currentunits = get(gca,'Units');
            set(gca, 'Units', 'Points');
            axpos = get(gca,'Position');
            set(gca, 'Units', currentunits);
            markerWidth = 0.01/diff(xlim)*axpos(3); % Calculate Marker width in points
            set(h, 'SizeData', markerWidth^2)
            line('XData', [0 aod_max_lim], 'YData', [0 aod_max_lim], 'LineStyle', '-','LineWidth', 1, 'Color','k')
            legend(varargin{2:end},'Location','northwest')
            xlabel('AERONET Measurement','FontSize',18)
            ylabel('AOD Retrieval','FontSize',18)
        end
        
        %legend({'Random Local Search','Coordinate Ascent','MCMC'},'Location','southeast')
        set(gca,'FontSize',18)
        export_fig(strcat('plots/',plot_name,'_',strjoin(varargin,'_')),'-png','-transparent','-r240')

    elseif strcmp(plot_name,'overlay')
        Date = varargin{1};
        Path = varargin{2};
        Orbit = varargin{3};
        Block = varargin{4};
        Location = varargin{5};
        Method = varargin{6};
        cmap = varargin{7};
        
        [aod_model, ~,~, ~, lon1,lat1] = load_aod(Date,Path,Orbit,Block,const,Method);
        [aod_aeronet, ~, ~, lon2,lat2] = load_aeronet(Date,Path,Block,Location,const);
        
        aod = [aod_model;aod_aeronet];
        %aod_min = min(aod);aod_max = max(aod);
        aod_min = 0.02;aod_max=0.16;
        
        cols = colormap(cmap);
        colsize = size(cols,1);
        map_model = round(1 + (aod_model - aod_min) / (aod_max-aod_min) .* (colsize-1));
        map_aeronet = round(1 + (aod_aeronet - aod_min) / (aod_max-aod_min) .* (colsize-1));
        
        map_model(map_model>256)=256;map_model(map_model<1)=1;
        map_aeronet(map_aeronet>256)=256;map_aeronet(map_aeronet<1)=1;
        
        h = scatter_patches(lon1,lat1,36,cols(map_model,:),'s','FaceAlpha',0.6,'EdgeColor','none');hold on
        scatter_patches(lon2,lat2,50,cols(map_aeronet,:),'o','EdgeColor',[0 0 0]);
        colorbar,caxis([aod_min aod_max])
        uistack(h,'bottom');
        
        plot_google_map('MapType','hybrid','ShowLabels',0,'Alpha',0.8)
        file_model = strcat('plots/',plot_name,'_',strjoin({Date,Location,Method},'_'));
        %export_fig(file_model,'-png','-transparent','-r240')
    
    elseif strcmp(plot_name,'resid')
        
        Date = varargin{1};
        Path = varargin{2};
        Orbit = varargin{3};
        Block = varargin{4};
        Method = varargin{5};
        iter = varargin{6};
        
        sample = load_cache(Date,Path,Orbit,Block,const,'sample',Method,0,0,1);
        
        resid = squeeze(sample.resid(:,:,iter))';
        resid(isinf(resid)) = NaN;
       
        C = nancov(resid);
        [R,~] = corrcov(C);
        imagesc(R),colorbar;
        set(gcf,'Position',[100, 100, 800, 600]);
        set(gca,'FontSize',18)
        export_fig(strcat('plots/',plot_name,'_',strjoin({Date,Method,num2str(iter)})),'-png','-transparent','-r240')
        
    elseif strcmp(plot_name,'stability')
       
        Location = varargin{1};
        test_delta = varargin{2};
        
        figure
        aod_max = 0;
        for p = 3:length(varargin)
            Method = varargin{p};

            [aoda,aodb] = load_aod_rep(Location,const,Method,test_delta);
            
            if test_delta == 0
                avga  = mean(aoda,2);
                stda = std(aoda,0,2); 
                errorbar(aodb,avga,-stda,stda,'MarkerEdgeColor',cols(p-2,:),'MarkerFaceColor',cols(p-2,:),'Color',cols(p-2,:),'LineStyle','none','Marker','x')
                aod_max = max([aod_max;avga+stda;aodb]); 
            else
                plot(repmat(aodb,1,6),aoda,'x')
                legend({'\Delta=0.001','\Delta=0.01','\Delta=0.05','\Delta=0.1','\Delta=1','\Delta=10'},'Location','southeast')
                aod_max = max([aod_max;aoda(:);aodb]);
                delta_all = [0.001,0.01,0.05,0.1,1,10];
                
                for r = 1:6
                    fprintf('delta is %f\n',delta_all(r))
                    fprintf('correlation with aeronet is %f\n',corr(aodb,aoda(:,r)))
                    fprintf('rms is %f\n',norm(aodb-aoda(:,r))/sqrt(length(aodb)))
                    fprintf('avg bias is %f\n',mean(aodb-aoda(:,r)))
                end
                
                
            end
            
            hold on
        end
        aod_max_lim = 1.05*aod_max;
        xlim([0,aod_max_lim]),ylim([0,1.05*aod_max_lim])
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
        export_fig(file_model,'-png','-transparent','-r240')
        
        
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
          
    end

end