function clean_product(Date)

    dir_radiance = fullfile('products/MI1B2T/',Date);
    if exist(dir_radiance,'dir')
        rmdir(dir_radiance,'s')
    else
        fprintf('%s not exists, not need to remove!\n',dir_radiance)
    end

end