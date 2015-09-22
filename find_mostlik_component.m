function [Component_Particle,Component_Num] = find_mostlik_component(reg,smart,x,y,ExtCroSect,CompSSA,kf,const,add_limit)

    cor1 = zeros(reg.num_reg_used,const.Band_Dim);
    cor2 = zeros(reg.num_reg_used,const.Band_Dim);
    ratio = zeros(reg.num_reg_used,const.Band_Dim);
    component = zeros(reg.num_reg_used,const.Band_Dim);
    const.Component_Particle = 1:21;
    const.Component_Num = length(const.Component_Particle);
    
    parfor p=1:reg.num_reg_used
        [regp,smartp] = extract_pixel(x(p),y(p),reg,smart,const,kf);
        [cor1(p,:),cor2(p,:),ratio(p,:),component(p,:)] = extract_cor2(regp,smartp,ExtCroSect,CompSSA,kf,const,add_limit);
    end
    
    sorted_table = sortrows(tabulate(component(:,1)),2);
    Component_Particle = int8(sorted_table(end-7:end,1));
    Component_Num = length(Component_Particle);
    
    disp(Component_Particle')

end