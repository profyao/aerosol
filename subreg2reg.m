function reg = subreg2reg(rad_1100,rad_275,Date,Path,Orbit,Block,r,const)

    if r<1100
        error('resolution not supported!')
    end
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    % Subregion used for HetSurf algorithm
    try 
        SubrUsed = hdfread(file_aerosol, 'SubregParamsAer', 'Fields', 'SubrUsed','Index',{[Block  1  1],[1  1  1],[1  const.XDim_r1100  const.YDim_r1100]}); %128x512
        AlgTypeFlag = hdfread(file_aerosol, 'RegParamsAer', 'Fields', 'AlgTypeFlag','Index',{[Block  1  1],[1  1  1],[1  const.XDim_r17600  const.YDim_r17600]}); %8x32
        NumAcceptSubr = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'NumAcceptSubr','Index',{[Block  1  1  1  1],[1  1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  const.Band_Dim  const.Cam_Dim]});
    catch ME
        rethrow(ME)
    end

    [I, J] = find(SubrUsed ~= 2);
    for ii = 1:length(I)
       rad_1100(I(ii), J(ii), :, :) = NaN;
       rad_275(I(ii)*4-3:4*I(ii),J(ii)*4-3:4*J(ii),:,:) = NaN;
    end
    
    XDim_r = const.XDim_r4400 * const.r4400/r;
    YDim_r = const.YDim_r4400 * const.r4400/r;

    RegSize = r/const.r1100;
    RegScale = const.r17600/r;
    
    reg.reg_is_used = false(XDim_r, YDim_r);
    reg.channel_is_used = false(XDim_r, YDim_r, const.Band_Dim, const.Cam_Dim);

    reg.mean_equ_ref = NaN * ones(XDim_r, YDim_r, const.Band_Dim, const.Cam_Dim);
    reg.min_equ_ref = NaN * ones(XDim_r, YDim_r, const.Band_Dim, const.Cam_Dim);
    
    if r > const.r1100
        reg.eof = NaN * ones(XDim_r, YDim_r,const.Band_Dim, const.Cam_Dim, const.Cam_Dim);
        reg.eigenvalue = zeros(XDim_r, YDim_r ,const.Band_Dim, const.Cam_Dim);
        reg.max_usable_eof = NaN * ones(XDim_r, YDim_r, const.Band_Dim);
    else
        reg.eof = NaN * ones(XDim_r, YDim_r,const.Cam_Dim, const.Cam_Dim);
        reg.eigenvalue = zeros(XDim_r, YDim_r , const.Cam_Dim);
        reg.max_usable_eof = NaN * ones(XDim_r, YDim_r);
    end
        
    %fprintf('start converting subreg radiance to reg!\n')

    for ii = 1:XDim_r
        
        for jj = 1:YDim_r
            
            ii_4400 = ceil(ii * r/const.r4400);
            jj_4400 = ceil(jj * r/const.r4400);
            
            subr_used_4400 = SubrUsed((const.RegSize*(ii_4400-1)+1):const.RegSize*ii_4400, (const.RegSize*(jj_4400-1)+1):const.RegSize*jj_4400)==2;            
            
            if  AlgTypeFlag(ceil(ii/RegScale), ceil(jj/RegScale)) == 3 && sum(subr_used_4400(:)) >= const.Config_min_het_subr_thresh
                
                if r == const.r1100
                    subrs_275 = rad_275(ii*4-3:ii*4,jj*4-3:jj*4,:);
                    subrs_275 = reshape(subrs_275, 4*4, const.Cam_Dim);
                end

                for band = 1:const.Band_Dim
                  
                    cam_is_used = reshape(NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)>0 ...
                        & NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)~=65533,1,const.Cam_Dim); % camera used
                    
                    subrs_1100 = rad_1100((RegSize*(ii-1)+1):RegSize*ii, (RegSize*(jj-1)+1):RegSize*jj, band, :);            
                    subrs_1100 = reshape(subrs_1100, RegSize*RegSize, const.Cam_Dim);
                    
                    % determine used cams                    
                    if r==1100
                        [sample_matrix,cam_is_used] = remove_nan(subrs_275,cam_is_used,const);
                    else
                        [sample_matrix,cam_is_used] = remove_nan(subrs_1100,cam_is_used,const);
                    end
                    
                    if isempty(sample_matrix)
                        continue
                    end
                    
                    num_cam_is_used = sum(cam_is_used);
                    
                    % Region averaged equivalent reflectances
                    reg.reg_is_used(ii, jj) = true;
                    reg.mean_equ_ref(ii, jj, band, :) = nanmean(subrs_1100, 1);
                    reg.mean_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                    reg.min_equ_ref(ii, jj, band, :) = min(subrs_1100, [], 1);
                    reg.min_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                    reg.channel_is_used(ii, jj, band, cam_is_used) = true;
                    
                    if r==1100 && band ~= const.Band_Red
                        continue
                    end
                    
                    % Calculate empirical orthogonal functions
                    sample_matrix = sample_matrix(:, cam_is_used);

                    [~,s,v] = svd(sample_matrix);
                    
                    d = diag(s).^2;
                    
                    if r > 1100
                        reg.eof(ii, jj, band, 1:num_cam_is_used, 1:num_cam_is_used) = v;
                        reg.eigenvalue(ii, jj, band, 1:length(d)) = d;
                    else
                        reg.eof(ii,jj,1:num_cam_is_used,1:num_cam_is_used) = v;
                        reg.eigenvalue(ii,jj,1:length(d)) = d;
                    end

                    cumd = cumsum(d);
                    max_usable = find(cumd >= const.Config_eigenvector_variance_thresh * cumd(end), 1, 'first');
                    
                    if r > 1100
                        reg.max_usable_eof(ii, jj, band) = max_usable;
                    else
                        reg.max_usable_eof(ii, jj) = max_usable;
                    end
                    
                end
                
            else
                %fprintf('.')
            end
        end
    end
    
    reg.ind_used = double(reg.reg_is_used);
    reg.num_reg_used = sum(reg.reg_is_used(:));
    reg.ind_used(reg.reg_is_used) = 1:reg.num_reg_used;
    
end
