function preprocess_batch(Location,const,r)

    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));
    
    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
        
    for i = 1:length(Dates)
        
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        %download product
        download_product(Date,Path,Orbit,const);
        
        %stage-1
        [rad_1100, ~, rad_275] = load_subreg(Date,Path,Orbit,Block,const);        
        %[rad_1100, rad_275] = load_cache(Date,Path,Orbit,Block,const,'rad_1100','rad_275');
        
        %stage-2
        reg = subreg2reg(rad_1100,rad_275,Date,Path,Orbit,Block,r,const);
        
        %stage-3
        smart = load_smart(Date,Path,Orbit,Block,const,0);
        %smart = load_cache(Date,Path,Orbit,Block,r,'smart');
        
        %save preprocessed data
        save2cache(Date,Path,Orbit,Block,r,reg,smart);
        
        %clean_product(Date);

    end

end