function [sample_matrix,cam_is_used] = remove_nan(subrs,cam_is_used,const)

    % Find reference cam for each band
    ref_cams = [5,4,6,3,7,2,8,1,9];
    for i = 1:const.Cam_Dim
        if ismember(ref_cams(i),find(cam_is_used))
            offset_equ_ref = subrs(:, ref_cams(i)); % references channel
            subreg_is_used = ~isnan(offset_equ_ref);
            num_subreg_used = sum(subreg_is_used);
            if num_subreg_used >= const.Config_min_het_subr_thresh
                break
            end
        else
            continue
        end
    end

    if num_subreg_used < const.Config_min_het_subr_thresh
        sample_matrix = [];
        cam_is_used = [];
    else
        [~, ndx] = min(offset_equ_ref);

        sample_matrix = subrs - repmat(subrs(ndx, :), 16, 1);
        sample_matrix = sample_matrix( ~ismember(1:16,ndx)' & subreg_is_used,:);

        cid = sum(isnan(sample_matrix),1)==0;
        cam_is_used = cam_is_used & cid;
    end

end