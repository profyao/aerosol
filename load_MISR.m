function [aod,frac,xid,yid,valid] = load_MISR(Date,Path,Orbit,Block,r,const)
        
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_F12_0022.hdf');
    tau0 = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepth', ...
    'Index',{[Block  1  1  const.Band_Green],[1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  1]});
    tau0 = double(tau0);
    
    RegScale = const.r17600/r;
    XDim_r = const.XDim_r4400 * const.r4400/r;
    YDim_r = const.YDim_r4400 * const.r4400/r;
    
    tau0_2d = kron(tau0, ones(RegScale)); 
    valid = tau0_2d ~= -9999;
    aod = tau0_2d(valid);
    [xid,yid] = ind2sub([XDim_r, YDim_r],find(valid));
    
    frac0 = hdfread(file_aerosol,'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepthFraction',...
        'Index',{[Block  1  1  const.Band_Green 1],[1 1 1 1 1],[1 const.XDim_r17600  const.YDim_r17600  1  5]});
    frac0 = double(squeeze(frac0));
    num_valid = sum(valid(:));
    frac = NaN*ones(num_valid,5);
    
    for i = 1:5
        tmp = kron(frac0(:,:,i),ones(RegScale));
        frac(:,i) = tmp(valid);
    end
    
    frac = frac';
        
    fprintf('MISR AOD from %s is loaded!\n',file_aerosol)
            
end