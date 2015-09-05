function [lon,lat] = get_coord(Path,Block,const)
    
    header_MIANCAGP_filename = const.header_MIANCAGP_filename;
    
    Path = num2str(Path,'%03d');
    
    dir_geo = fullfile('products/MIANCAGP'); 
    file_geo = strcat(dir_geo,'/',header_MIANCAGP_filename,Path,'_F01_24.hdf');
   
    lat = hdfread(file_geo, 'Standard', 'Fields', 'GeoLatitude', 'Index', {[Block 1  1],[1  1  1],[1  const.XDim_r1100  const.YDim_r1100]});
    lon = hdfread(file_geo, 'Standard', 'Fields', 'GeoLongitude', 'Index', {[Block 1  1],[1  1  1],[1  const.XDim_r1100  const.YDim_r1100]});
    
end