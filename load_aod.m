function [aod,theta,xid,yid,lon,lat] = load_aod(Date,Path,Orbit,Block,r,const,Opt)

    RegSize = r/const.r1100;
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_F12_0022.hdf');
    ExtCroSect = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
    'Spectral extinction cross section', 'FirstRecord',1 ,'NumRecords', const.Model_ComponentDim);
    ExtCroSect = ExtCroSect{1}; % RH and band dependent

    [lon,lat] = get_coord(Path,Block,const);
    [lon,lat] = conv_coord(lon,lat,RegSize,const);

    if strcmp(Opt,'MISR')

        [aod,theta,xid,yid,valid] = load_MISR(Date,Path,Orbit,Block,r,const);
        lon = lon(valid);lat=lat(valid);
        return

    %elseif strcmp(Opt,'CD-random') || strcmp(Opt,'CD') || strcmp(Opt,'MCMC')  || strcmp(Opt,'CD-random-noprior') || strcmp(Opt,'CD-noprior')
    else
        [reg,sample] = load_cache(Date,Path,Orbit,Block,r,'reg','sample',Opt);
        [xid,yid] = find(reg.reg_is_used);
    %else
    %    error('no aod data source is specified!\n')
    end

    if strcmp(Opt,'MCMC')
        iter = size(sample.tau,2);
        tau = mean(sample.tau(:,200:50:end),2);
        theta = squeeze(mean(sample.theta(:,:,200:50:end),3));
        %theta = squeeze(sample.theta(:,:,1));
    else
        tau = sample.tau;
        theta = sample.theta;
    end

    aod = extract_aod_4band(tau,theta,reg.num_reg_used,ExtCroSect,const);

    lon = lon(reg.reg_is_used);
    lat = lat(reg.reg_is_used);
        
   

    
end
