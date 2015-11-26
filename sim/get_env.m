function [ExtCroSect, CompSSA, smart, reg] = get_env(Date,Path,Orbit,Block,r,const)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');
    
    ExtCroSect = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
            'Spectral extinction cross section', 'FirstRecord',1 ,'NumRecords', const.Model_ComponentDim);
    ExtCroSect = ExtCroSect{1}; % RH and band dependent
    CompSSA = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
        'Spectral single scattering albedo', 'FirstRecord',1 ,'NumRecords', const.Model_ComponentDim);
    CompSSA = CompSSA{1};
    
    [smart,reg] = load_cache(Date,Path,Orbit,Block,r,'smart','reg');


end