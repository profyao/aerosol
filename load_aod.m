function [aod,xid,yid,lon,lat] = load_aod(Date,Path,Orbit,Block,const,Opt)
    
    try
        [lon,lat] = get_coord(Path,Block,const);
        [lon,lat] = conv_coord(lon,lat,const.RegSize,const);
            
        if strcmp(Opt,'MISR')
            
            dir_aerosol = fullfile('products/MIL2ASAE/',Date);
            file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_F12_0022.hdf');
            tau0 = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepth', ...
        'Index',{[Block  1  1  const.Band_Green],[1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  1]});
            tau0 = double(tau0);
            tau0_2d = kron(tau0, ones(const.RegScale)); 
            valid = tau0_2d ~= -9999;
            aod = tau0_2d(valid);
            [xid,yid] = ind2sub([const.XDim_r, const.YDim_r],find(valid));
            lon = lon(valid);
            lat=lat(valid);
            fprintf('MISR AOD from %s is loaded!\n',file_aerosol)
            return

        elseif strcmp(Opt,'CD-random') || strcmp(Opt,'MCMC') || strcmp(Opt,'CD')

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
        else
            aod = sample.tau(:,end);
        end

        lon = lon(reg.reg_is_used);
        lat = lat(reg.reg_is_used);
        
    catch
        aod = [];
        xid = [];
        yid =[];
        lon = [];
        lat = [];
    end
    
end
