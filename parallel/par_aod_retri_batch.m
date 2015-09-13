function par_aod_retri_batch(Method,Location,const)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    N = length(Dates);
    
    for i = 8:8
        
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        [sample,error_flag] = par_aod_retri(Date,Path,Orbit,Block,Method,const);
        
        if (error_flag == 1)
            fprintf(strcat(Date,'_P',num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_B',num2str(Block),' is skipped!\n'))
            continue
        else
            save2cache(Date,Path,Orbit,Block,const,sample,Method)
        end

    end

end