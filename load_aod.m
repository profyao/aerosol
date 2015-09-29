function [aod,theta,xid,yid,lon,lat] = load_aod(Date,Path,Orbit,Block,const,Opt)
    
    try
        [lon,lat] = get_coord(Path,Block,const);
        [lon,lat] = conv_coord(lon,lat,const.RegSize,const);
            
        if strcmp(Opt,'MISR')
            
            [aod,theta,xid,yid,valid] = load_MISR(Date,Path,Orbit,Block,const);
            lon = lon(valid);lat=lat(valid);
            return

        elseif strcmp(Opt,'CD-random') || strcmp(Opt,'MCMC') || strcmp(Opt,'CD') || strcmp(Opt,'CD-random-noprior')

            [reg,sample] = load_cache(Date,Path,Orbit,Block,const,'reg','sample',Opt,0,0,1); %Method,kf,dy,par
            
        elseif strcmp(Opt,'CD-random-par')
            
            [reg,sample] = load_cache(Date,Path,Orbit,Block,const,'reg','sample','CD-random',0,0,1); %Method,kf,dy,par
            
        elseif strcmp(Opt,'CD-random-kf')
            
            [reg,sample] = load_cache(Date,Path,Orbit,Block,const,'reg','sample','CD-random',1,0,1);
            
        elseif strcmp(Opt,'CD-random-dy')
            
            [reg,sample] = load_cache(Date,Path,Orbit,Block,const,'reg','sample','CD-random',0,1,1);
            
        elseif strcmp(Opt,'CD-random-nopar')
            
            [reg,sample] = load_cache(Date,Path,Orbit,Block,const,'reg','sample','CD-random',0,0,0);
            
        else
            error('no aod data source is specified!\n')
        end
        
        [xid,yid] = find(reg.reg_is_used);
        
        if strcmp(Opt,'MCMC')
            iter = size(sample.tau,2);
            aod = mean(sample.tau(:,ceil(iter/2):end),2);
            theta = squeeze(mean(sample.theta(:,:,ceil(iter/2):end),3));
        else
            aod = sample.tau;
            theta = sample.theta;
        end

        lon = lon(reg.reg_is_used);
        lat = lat(reg.reg_is_used);
        
        %varargout{1} = reg;
        %varargout{2} = sample;
        
    catch
        aod = [];
        theta = [];
        xid = [];
        yid =[];
        lon = [];
        lat = [];
        %varargout={};
    end
    
end
