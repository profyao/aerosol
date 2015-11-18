function dat_n = add_noise(dat,k)
    
    dim = size(dat);
    
    tmp = reshape(dat, [dim(1)*dim(2),dim(3),dim(4)]);
    dat_avg = squeeze(nanmean(tmp));
    dat_avg = dat_avg(:);
    
    noise = k * randn(dim(1)*dim(2), dim(3)*dim(4)) .* repmat(dat_avg,[1,dim(1)*dim(2)])';
    
    noise = reshape(noise,[dim(1), dim(2), dim(3), dim(4)]);
    
    dat_n = dat + noise;
    
    dat_n(dat_n<=0) = 0.01;

end