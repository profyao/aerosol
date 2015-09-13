function reg = subreg2reg2(subreg,Date,Path,Orbit,Block,const)
    
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

    reg.eof = NaN * ones(XDim_r, YDim_r,Cam_Dim, Cam_Dim);
    reg.eigenvalue = NaN * ones(XDim_r, YDim_r ,Cam_Dim);
    reg.max_usable_eof = NaN * ones(XDim_r, YDim_r);
    
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
                num_subreg_used = sum(subreg_used);
                [~, ndx] = nanmin(offset_equ_ref); % the darkest subregion, the alogrithm seems to be too sensitive to the choice of the reference subregion
                cam_is_used_allband = true(1,Cam_Dim);
                sample_matrix_allband = NaN*ones(num_subreg_used*Band_Dim,Cam_Dim);

                for band = 1:Band_Dim
                    
                    cam_is_used = reshape(NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)~=0 ...
                        & NumAcceptSubr(ceil(ii/RegScale), ceil(jj/RegScale), band, :)~=65533,Cam_Dim,1); % camera used
                    num_cam_is_used = sum(cam_is_used);
                    
                    cam_is_used_allband = cam_is_used .* cam_is_used_allband;
                    
                    if num_cam_is_used>0
                        
                        reg.channel_is_used(ii, jj, band, cam_is_used) = 1;
                        reg.num_cam_used(ii, jj, band) = num_cam_is_used;
                        reg.mean_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                        reg.num_subreg(ii, jj, band, ~cam_is_used) = 0;
                        reg.min_equ_ref(ii, jj, band, ~cam_is_used) = NaN;

                        sample_matrix = reshape(subrs(subreg_used, band, :),num_subreg_used,Cam_Dim);
                        reduced_sample_matrix = sample_matrix - repmat(reshape(subrs(ndx, band, :),1,Cam_Dim), num_subreg_used, 1);               

                        sample_matrix_allband((band-1)*num_subreg_used+1:band*num_subreg_used,:) = reduced_sample_matrix;
                   
                    else
                        fprintf('%d,%d: no camera meets retrieval quality!\n',ii,jj)
                    end
                    
                end
                
                sample_matrix_allband = sample_matrix_allband(sample_matrix_allband(:,1)~=0,cam_is_used_allband);
                num_cam_is_used_allband = sum(cam_is_used_allband);

                [~,s,v] = svd(sample_matrix_allband);
                d = diag(s).^2;
                
                reg.eof(ii, jj, 1:num_cam_is_used_allband,1:num_cam_is_used_allband) = v;
                reg.eigenvalue(ii, jj, 1:num_cam_is_used_allband) = d;

                cumd = cumsum(d);
                reg.max_usable_eof(ii, jj) = max(Config_first_eigenvalue_for_eofs, ...
                    find(cumd >= Config_eigenvector_variance_thresh * cumd(end), 1, 'first'));
                        
            else
                %fprintf('dark water or not heterogeneous algorithm!\n')
            end
        end
    end
    
    reg.ind_used = double(reg.reg_is_used);
    reg.num_reg_used = sum(reg.reg_is_used(:));
    reg.ind_used(reg.reg_is_used) = 1:reg.num_reg_used;

end
