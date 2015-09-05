function reg = subreg2reg(subreg,Date,Path,Orbit,Block,const)
    
    header_MIL2ASAE_filename = const.header_MIL2ASAE_filename;
    XDim_r1100 = const.XDim_r1100;
    YDim_r1100 = const.YDim_r1100;
    XDim_r = const.XDim_r;
    YDim_r = const.YDim_r;
    XDim_r17600 = const.XDim_r17600;
    YDim_r17600 = const.YDim_r17600;
    Band_Dim = const.Band_Dim;
    Band_Green = const.Band_Green;
    Cam_AN = const.Cam_AN;
    Cam_Dim = const.Cam_Dim;
    RegScale = const.RegScale;
    RegSize = const.RegSize;
    Config_min_het_subr_thresh = const.Config_min_het_subr_thresh;
    Config_first_eigenvalue_for_eofs = const.Config_first_eigenvalue_for_eofs;
    Config_eigenvector_variance_thresh = const.Config_eigenvector_variance_thresh;

    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    % Subregion used for HetSurf algorithm
    try 
        SubrUsed = hdfread(file_aerosol, 'SubregParamsAer', 'Fields', 'SubrUsed','Index',{[Block  1  1],[1  1  1],[1  XDim_r1100  YDim_r1100]}); %128x512
        NumAcceptSubr = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'NumAcceptSubr','Index',{[Block  1  1  1  1],[1  1  1  1  1],[1  XDim_r17600  YDim_r17600  Band_Dim  Cam_Dim]});
        AlgTypeFlag = hdfread(file_aerosol, 'RegParamsAer', 'Fields', 'AlgTypeFlag','Index',{[Block  1  1],[1  1  1],[1  XDim_r17600  YDim_r17600]}); %8x32, 17.6km
    catch ME
        rethrow(ME)
    end

    [I, J] = find(SubrUsed ~= 2);
    for ii = 1:length(I)
       subreg(I(ii), J(ii), :, :) = NaN; 
    end
    
    reg.reg_is_used = false(XDim_r, YDim_r);
    reg.channel_is_used = false(XDim_r, YDim_r, Band_Dim, Cam_Dim);
    reg.num_cam_used =  zeros(XDim_r, YDim_r, Band_Dim);

    reg.mean_equ_ref = NaN * ones(XDim_r, YDim_r, Band_Dim, Cam_Dim);
    reg.min_equ_ref = NaN * ones(XDim_r, YDim_r, Band_Dim, Cam_Dim);
    reg.num_subreg = NaN * ones(XDim_r, YDim_r, Band_Dim, Cam_Dim);

    reg.eof = NaN * ones(XDim_r, YDim_r,Band_Dim, Cam_Dim, Cam_Dim);
    reg.eigenvalue = NaN * ones(XDim_r, YDim_r ,Band_Dim, Cam_Dim);
    reg.max_usable_eof = NaN * ones(XDim_r, YDim_r, Band_Dim);
    
    fprintf('convert subreg to reg!\n')

    for ii = 1:XDim_r
        
        for jj = 1:YDim_r
            
            subr_used = SubrUsed((RegSize*(ii-1)+1):RegSize*ii, (RegSize*(jj-1)+1):RegSize*jj)==2;        
            
            if  AlgTypeFlag(ceil(ii/RegScale), ceil(jj/RegScale)) == 3 && sum(subr_used(:)) >= Config_min_het_subr_thresh
                
                %disp([ii,jj]);
                subrs = subreg((RegSize*(ii-1)+1):RegSize*ii, (RegSize*(jj-1)+1):RegSize*jj, :, :);            
                subrs = reshape(subrs, RegSize^2, Band_Dim, Cam_Dim);
                
                % Region averaged equivalent reflectances
                reg.reg_is_used(ii, jj) = true;
                reg.mean_equ_ref(ii, jj, :, :) = squeeze(nanmean(subrs, 1));            
                reg.num_subreg(ii, jj, :, :) = squeeze(sum(~isnan(subrs), 1));
                reg.min_equ_ref(ii, jj, :, :) = squeeze(nanmin(subrs, [], 1));

                % Calculate empirical orthogonal functions
                offset_equ_ref = subrs(:, Band_Green, Cam_AN); % references channel
                subreg_used = ~isnan(offset_equ_ref);
                [~, ndx] = nanmin(offset_equ_ref); % the darkest subregion, the alogrithm seems to be too sensitive to the choice of the reference subregion

                for band = 1:Band_Dim                
                    cam_is_used = reshape(NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)~=0 ...
                        & NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)~=65533,Cam_Dim,1); % camera used
                    
                    if sum(cam_is_used)>0
                        reg.channel_is_used(ii, jj, band, cam_is_used) = 1;
                        reg.num_cam_used(ii, jj, band) = sum(cam_is_used);
                        reg.mean_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                        reg.num_subreg(ii, jj, band, ~cam_is_used) = 0;
                        reg.min_equ_ref(ii, jj, band, ~cam_is_used) = NaN;

                        sample_matrix = reshape(subrs(subreg_used, band, cam_is_used),sum(subreg_used),sum(cam_is_used));

                        reduced_sample_matrix = sample_matrix - repmat(squeeze(subrs(ndx, band, cam_is_used))', size(sample_matrix, 1), 1);               

                        scatter_matrix = reduced_sample_matrix'*reduced_sample_matrix;

                        [v,d] = eig(scatter_matrix); 
                        % eigenvalue decomposition: scatter_matrix*v = v*d
                        [d, IX] = sort(diag(d), 'descend');
                        v = v(:, IX);

                        reg.eof(ii, jj, band, 1:reg.num_cam_used(ii,jj,band), 1:reg.num_cam_used(ii,jj,band)) = v;
                        reg.eigenvalue(ii, jj, band, 1:reg.num_cam_used(ii,jj,band)) = d;

                        cumd = cumsum(d);
                        reg.max_usable_eof(ii, jj, band) = max(Config_first_eigenvalue_for_eofs, ...
                            find(cumd >= Config_eigenvector_variance_thresh * cumd(end), 1, 'first'));
                    else
                        fprintf('%d,%d: no camera meets retrieval quality!\n',ii,jj)
                    end
                    
                end
            else
                %fprintf('dark water or not heterogeneous algorithm!\n')
            end
        end
    end
    
    reg.ind_used = double(reg.reg_is_used);
    reg.num_reg_used = sum(reg.reg_is_used(:));
    reg.ind_used(reg.reg_is_used) = 1:reg.num_reg_used;

end
