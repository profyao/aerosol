function download_product(Date,Path,Orbit,const)

    Cam_Name = const.Cam_Name;
    header_MI1B2T_url = const.header_MI1B2T_url;
    header_MI1B2T_filename = const.header_MI1B2T_filename;
    header_MIL2ASAE_url = const.header_MIL2ASAE_url;
    header_MIL2ASAE_filename = const.header_MIL2ASAE_filename;
    Cam_Dim = const.Cam_Dim;
    
    Cam_Name0 = Cam_Name;
    header_MI1B2T_url0 = header_MI1B2T_url;
    header_MI1B2T_filename0 = header_MI1B2T_filename;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    % MI1B2T download
    disp('MI1B2T..');
    
    dir_radiance = fullfile('products/MI1B2T/',Date);
    if ~exist(dir_radiance,'dir')
        mkdir(dir_radiance)
    else
        fprintf('%s exists, continue downloading...\n',dir_radiance)
    end
  
    parfor i = 1:Cam_Dim
        cam = Cam_Name0{i};
        file_name = strcat(header_MI1B2T_filename0,Path,'_O',Orbit,'_',cam,'_F03_0024.hdf');
        file_url = strcat(header_MI1B2T_url0, Date, '/', file_name);
        local_file_name = strcat(dir_radiance,'/', file_name);
        if ~exist(local_file_name,'file')
            fprintf('download from %s to %s ...\n',file_url,local_file_name)
            urlwrite(file_url, local_file_name);
        else
            fprintf('%s already exists, continue to next cam...\n',local_file_name)
            continue
        end
    end
    
    disp('MIL2ASAE..');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    if ~exist(dir_aerosol,'dir')
        mkdir(dir_aerosol)
    else
        fprintf('%s exists, continue downloading...\n',dir_aerosol)
    end
    
    file_name = strcat(header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf'); 
    file_url = strcat(header_MIL2ASAE_url,Date,'/',file_name);
    local_file_name = strcat(dir_aerosol,'/', file_name);
    
    if ~exist(local_file_name,'file')
        fprintf('download from %s to %s ...\n',file_url,local_file_name)
        urlwrite(file_url, local_file_name);
    else
        sprintf('%s already exists!\n',local_file_name);
    end
    
    disp('MIANCAGP..');
    
    dir_geo = fullfile('products/MIANCAGP');
    if ~exist(dir_geo,'dir')
        mkdir(dir_geo)
    else
        fprintf('%s exists, continue downloading...\n',dir_geo)
    end
    
    file_name = strcat(const.header_MIANCAGP_filename,Path,'_F01_24.hdf');
    file_url = {};
    file_url{1} = strcat(const.header_MIANCAGP_url1,file_name);
    file_url{2} = strcat(const.header_MIANCAGP_url2,file_name);
    local_file_name = strcat(dir_geo,'/', file_name);
    attempt = 1;
    
    if ~exist(local_file_name,'file')
        while attempt <= 2
            try
                urlwrite(file_url{attempt}, local_file_name);
                fprintf('download from %s to %s ...\n',file_url{attempt},local_file_name)
                break
            catch
                disp('An error occurred while retrieving information from the internet...')
            end
            attempt = attempt + 1;
        end     
    else
        sprintf('%s already exists!\n',local_file_name);
    end

end