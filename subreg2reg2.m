function reg = subreg2reg2(subreg,Date,Path,Orbit,Block,const)
    
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

    reg.eof_allband = NaN * ones(const.XDim_r, const.YDim_r,const.Cam_Dim, const.Cam_Dim);
    reg.eigenvalue_allband = zeros(const.XDim_r, const.YDim_r ,const.Cam_Dim);
    reg.max_usable_eof_allband = NaN * ones(const.XDim_r, const.YDim_r);
        
    fprintf('convert subreg to reg!\n')

    for ii = 1:const.XDim_r
        
        for jj = 1:const.YDim_r
            
            subr_used = SubrUsed((const.RegSize*(ii-1)+1):const.RegSize*ii, (const.RegSize*(jj-1)+1):const.RegSize*jj)==2;            
            
            if  AlgTypeFlag(ceil(ii/const.RegScale), ceil(jj/const.RegScale)) == 3 && sum(subr_used(:)) >= const.Config_min_het_subr_thresh
                
                cam_is_used_allband = true(1,const.Cam_Dim);
                sample_matrix_allband = NaN*ones(const.sample_size*const.Band_Dim,const.Cam_Dim);

                for band = 1:const.Band_Dim
                  
                    cam_is_used = reshape(NumAcceptSubr(ceil(ii/const.RegScale), ceil(jj/const.RegScale), band, :)>0 ...
                        & NumAcceptSubr(ceil(ii/const.RegScale), ceil(jj/const.RegScale), band, :)~=65533,1,const.Cam_Dim); % camera used

                    num_cam_is_used = sum(cam_is_used);
                    
                    if num_cam_is_used>=const.min_cam_used
                        
                        cam_is_used_allband = cam_is_used & cam_is_used_allband;
                                               
                        subrs = subreg((const.RegSize*(ii-1)+1):const.RegSize*ii, (const.RegSize*(jj-1)+1):const.RegSize*jj, band, :);            
                        subrs = reshape(subrs, const.sample_size, const.Cam_Dim);
                
                        % Find reference cam
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
                        
                        %if i~=1
                        %    fprintf('%d,%d in band %d: ref cam is %d\n',ii,jj,band,ref_cams(i))                           
                        %end
                
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
                        reg.min_equ_ref(ii, jj, band, :) = min(subrs, [], 1);
                        reg.min_equ_ref(ii, jj, band, ~cam_is_used) = NaN;
                        reg.num_subreg_used(ii, jj, band, :) = sum(~isnan(subrs));
                        reg.num_subreg_used(ii, jj, band, ~cam_is_used) = 0;
                        reg.channel_is_used(ii, jj, band, cam_is_used) = true;
                        reg.num_cam_used(ii, jj, band) = num_cam_is_used;

                        reduced_sample_matrix = subrs - repmat(subrs(ndx, :),const.sample_size,1);                        
                        sample_matrix_allband((band-1)*const.sample_size+1:band*const.sample_size,:) = reduced_sample_matrix;
                        
                        % per band EOF
                        
                        sample_matrix = subrs(subreg_used_ref, cam_is_used);
                        reduced_sample_matrix = sample_matrix - repmat(subrs(ndx, cam_is_used),num_subreg_used_ref, 1);                                       
                        reduced_sample_matrix = reduced_sample_matrix(reduced_sample_matrix(:,1)~=0,:);
                        
                        [~,s,v] = svd(reduced_sample_matrix);
                        d = diag(s).^2;

                        reg.eof(ii, jj, band, 1:num_cam_is_used, 1:num_cam_is_used) = v;
                        reg.eigenvalue(ii, jj, band, 1:length(d)) = d;

                        cumd = cumsum(d);
                        max_usable = find(cumd >= const.Config_eigenvector_variance_thresh * cumd(end), 1, 'first');
                        if isempty(max_usable)
                            reg.max_usable_eof(ii, jj, band) = 1;
                            fprintf('%d,%d, band %d: eof is incorrect!\n',ii,jj,band)
                        else
                            reg.max_usable_eof(ii, jj, band) = max_usable;
                        end
                        
                        
                    else
                        fprintf('%d,%d, band %d: fewer than two cams available!\n',ii,jj,band)
                    end
                    
                end
                
                rid = sample_matrix_allband(:,1)~=0 & ~isnan(sample_matrix_allband(:,ref_cams(i)));
                sample_matrix_allband = sample_matrix_allband(rid,cam_is_used_allband);
                num_cam_is_used_allband = sum(cam_is_used_allband);
                
                [~,s,v] = svd(sample_matrix_allband);

                d = diag(s).^2;

                reg.eof_allband(ii, jj, 1:num_cam_is_used_allband,1:num_cam_is_used_allband) = v;
                reg.eigenvalue_allband(ii, jj, 1:length(d)) = d;
                cumd = cumsum(d);                
                max_usable = find(cumd >= const.Config_eigenvector_variance_thresh * cumd(end), 1, 'first');
                if isempty(max_usable)
                    reg.max_usable_eof_allband(ii, jj) = 1;
                    fprintf('%d,%d, all band eof is incorrect!\n',ii,jj)
                else
                    reg.max_usable_eof_allband(ii, jj) = max_usable;
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
