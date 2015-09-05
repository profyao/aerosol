function preprocess_batch(Location,const)

    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));
    
    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    N = length(Dates);
    
    for i = 1:N
        
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        download_product(Date,Path,Orbit,const);
        subreg = load_subreg(Date,Path,Orbit,Block,const);
        reg = subreg2reg(subreg,Date,Path,Orbit,Block,const);
        smart = load_smart(Date,Path,Orbit,Block,const);
        save2cache(Date,Path,Orbit,Block,const,subreg,reg,smart);
        clean_product(Date);

    end

end