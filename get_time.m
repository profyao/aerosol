function time = get_time(Date,Path,Orbit,Block,const)

    header_MIL2ASAE_filename = const.header_MIL2ASAE_filename;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    if exist(file_aerosol,'file')
        time = hdfread(file_aerosol,'PerBlockMetadataTime','Fields','BlockCenterTime','FirstRecord',1,'NumRecords',139);
        time = time{1};
        time = time(12:19,Block)';
    else
        fprintf('%s not exist!\n',file_aerosol)
    end
    
end
