function preprocess_batch(Location,const,add_limit)

    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));
    
    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    %N = length(Dates);
    
    parfor i = [8,9,11,12,14]
        
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        %download product
        %download_product(Date,Path,Orbit,const);
        
        %stage-1
        %subreg = load_subreg(Date,Path,Orbit,Block,const);        
        %subreg = load_cache(Date,Path,Orbit,Block,const,'subreg');
        
        %stage-2
        %reg = subreg2reg2(subreg,Date,Path,Orbit,Block,const);
        
        %stage-3
        smart = load_smart(Date,Path,Orbit,Block,const,add_limit);
        
        %save preprocessed data
        %save2cache(Date,Path,Orbit,Block,const,subreg,reg,smart);
        save2cache(Date,Path,Orbit,Block,const,smart);
        
        %clean_product(Date);

    end

end