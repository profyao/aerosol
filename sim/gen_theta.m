function theta = gen_theta(alpha, XDim_r, YDim_r)
    
    siz = length(alpha);
    len = XDim_r*YDim_r;
    theta = zeros(siz,len);
    tmp = zeros(siz,1);
    
    for i = 1:len
        for c = 1:siz
            tmp(c) = gamrnd(alpha(c),1);
        end
        theta(:,i) = tmp/sum(tmp);
    end

end