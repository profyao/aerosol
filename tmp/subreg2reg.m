function reg = subreg2reg(subreg,Date,Path,Orbit,Block,const)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    % Subregion used for HetSurf algorithm
    try 
        SubrUsed = hdfread(file_aerosol, 'SubregParamsAer', 'Fields', 'SubrUsed','Index',{[Block  1  1],[1  1  1],[1  const.XDim_r1100  const.YDim_r1100]}); %128x512
        NumAcceptSubr = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'NumAcceptSubr','Index',{[Block  1  1  1  1],[1  1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  const.Band_Dim  const.Cam_Dim]});
        AlgTypeFlag = hdfread(file_aerosol, 'RegParamsAer', 'Fields', 'AlgTypeFlag','Index',{[Block  1  1],[1  1  1],[1  const.XDim_r17600  const.YDim_r17600]}); %8x32, 17.6km
    catch ME
        rethrow(ME)
    end

    [I, J] = find(SubrUsed ~= 2);
    for ii = 1:length(I)
       subreg(I(ii), J(ii), :, :) = NaN; 
    end
    
    reg.reg_is_used = false(const.XDim_r, const.YDim_r);
    reg.channel_is_used = false(const.XDim_r, const.YDim_r, const.Band_Dim, const.Cam_Dim);
    reg.num_cam_used =  zeros(const.XDim_r, const.YDim_r, const.Band_Dim);

    reg.mean_equ_ref = NaN * ones(const.XDim_r, const.YDim_r, const.Band_Dim, const.Cam_Dim);
    reg.min_equ_ref = NaN * ones(const.XDim_r, const.YDim_r, const.Band_Dim, const.Cam_Dim);
    reg.num_subreg_used = NaN * ones(const.XDim_r, const.YDim_r, const.Band_Dim, const.Cam_Dim);

    reg.eof = NaN * ones(const.XDim_r, const.YDim_r,const.Band_Dim, const.Cam_Dim, const.Cam_Dim);
    reg.eigenvalue = zeros(const.XDim_r, const.YDim_r ,const.Band_Dim, const.Cam_Dim);
    reg.max_usable_eof = NaN * ones(const.XDim_r, const.YDim_r, const.Band_Dim);
    
    fprintf('prepare to convert subreg to reg!\n')

    for ii = 1:const.XDim_r
        
        for jj = 1:const.YDim_r
            
            subr_used = SubrUsed((const.RegSize*(ii-1)+1):const.RegSize*ii, (const.RegSize*(jj-1)+1):const.RegSize*jj)==2;        
            
            if  AlgTypeFlag(ceil(ii/const.RegScale), ceil(jj/const.RegScale)) == 3 && sum(subr_used(:)) >= const.Config_min_het_subr_thresh
                
                for band = 1:const.Band_Dim
                    
                    cam_is_used = reshape(NumAcceptSubr(ceil(ii/const.RegScale), ceil(jj/const.RegScale), band, :)>0 ...
                        & NumAcceptSubr(ceil(ii/const.RegScale), ceil(jj/const.RegScale), band, :)~=65533,1,const.Cam_Dim); % camera used

                    num_cam_is_used = sum(cam_is_used);
                    
                    if num_cam_is_used>=const.min_cam_used
                        
                        subrs = subreg((const.RegSize*(ii-1)+1):const.RegSize*ii, (const.RegSize*(jj-1)+1):const.RegSize*jj, band, :);            
                        subrs = reshape(subrs, const.sample_size, const.Cam_Dim);
                        
                        ref_cams = [5,4,6,3,7,2,8,1,9];
                        for i = 1:const.Cam_Dim
                            if ismember(ref_cams(i),find(cam_is_used))
                                offset_equ_ref = subrs(:, ref_cams(i)); % references channel
                                subreg_used_ref = ~isnan(offset_equ_ref);
                                num_subreg_used_ref = sum(subreg_used_ref);
                                if num_subreg_used_ref >= const.Config_min_het_subr_thresh
                                    break
                                end
                            else
                                continue
                            end
                        end
                        
                        if num_subreg_used_ref < const.Config_min_het_subr_thresh
                            fprintf('%d,%d, band %d: has fewer than 4 subregs that are not NaN, continue to next region!\n',ii,jj,band)
                            continue
                        end

                        % Calculate empirical orthogonal functions                    
                        [~, ndx] = min(offset_equ_ref);
                        %fprintf('%d,%d in band %d: darkest subreg is %d\n',ii,jj,band,ndx)
            
                        % Region averaged equivalent reflectances
                        reg.reg_is_used(ii, jj) = true;
                        reg.mean_equ_ref(ii, jj, band, :) = nanmean(subrs, 1);
                        reg.mean_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                        reg.min_equ_ref(ii, jj, band, :) = nanmin(subrs, [], 1);
                        reg.min_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                        reg.num_subreg_used(ii, jj, band, :) = squeeze(sum(~isnan(subrs), 1));
                        reg.num_subreg_used(ii, jj, band, ~cam_is_used) = 0;
                        reg.channel_is_used(ii, jj, band, cam_is_used) = true;
                        reg.num_cam_used(ii, jj, band) = num_cam_is_used;

               
                        % Calculate empirical orthogonal functions
                        sample_matrix = subrs(subreg_used_ref, cam_is_used);

                        reduced_sample_matrix = sample_matrix - repmat(subrs(ndx, cam_is_used),num_subreg_used_ref, 1);               

                        scatter_matrix = reduced_sample_matrix'*reduced_sample_matrix;

                        [v,d] = eig(scatter_matrix); 
                        % eigenvalue decomposition: scatter_matrix*v = v*d
                        [d, IX] = sort(diag(d), 'descend');
                        v = v(:, IX);

                        reg.eof(ii, jj, band, 1:num_cam_is_used, 1:num_cam_is_used) = v;
                        reg.eigenvalue(ii, jj, band, 1:num_cam_is_used) = d;

                        cumd = cumsum(d);
                        reg.max_usable_eof(ii, jj, band) = max(const.Config_first_eigenvalue_for_eofs, ...
                            find(cumd >= const.Config_eigenvector_variance_thresh * cumd(end), 1, 'first'));
                        
                    else
                        fprintf('%d,%d, band %d: fewer than two cams available!\n',ii,jj,band)
                    end
                    
                end
            else
                %fprintf('dark water or not heterogeneous algorithm!\n')
            end
        end
    end
    
    reg.num_reg_used = sum(reg.reg_is_used(:));

end
