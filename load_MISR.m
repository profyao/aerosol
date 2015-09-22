function [aod,theta,xid,yid,valid] = load_MISR(Date,Path,Orbit,Block,const)
        
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_F12_0022.hdf');
    tau0 = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepth', ...
    'Index',{[Block  1  1  const.Band_Green],[1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  1]});
    tau0 = double(tau0);
    tau0_2d = kron(tau0, ones(const.RegScale)); 
    valid = tau0_2d ~= -9999;
    aod = tau0_2d(valid);
    [xid,yid] = ind2sub([const.XDim_r, const.YDim_r],find(valid));
    theta = [];
    fprintf('MISR AOD from %s is loaded!\n',file_aerosol)
            
end