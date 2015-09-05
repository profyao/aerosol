function [cor1,cor2,ratio,component] = find_mostlik_component(Date,Path,Orbit,Block,band_ind,const)

    [reg,smart] = load_cache(Date,Path,Orbit,Block,const,'reg','smart');
    [x,y] = find(reg.reg_is_used);
    cor1 = NaN*ones(reg.num_reg_used,1);
    cor2 = NaN*ones(reg.num_reg_used,1);
    ratio = NaN*ones(reg.num_reg_used,1);
    component = NaN*ones(reg.num_reg_used,1);
    
    parfor p=1:reg.num_reg_used
        [corp,cor2(p),ratiop] = extract_cor(reg,smart,'SS',x(p),y(p),band_ind,const);
        [tmp,id1] = max(corp);
        [cor1(p),id2] = max(tmp);
        component(p) = id2;
        ratio(p) = ratiop(id1(id2),id2);
    end

end