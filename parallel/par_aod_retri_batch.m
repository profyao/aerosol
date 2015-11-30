function par_aod_retri_batch(Method,Location,r,core,const)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    N = length(Dates);
    
    %for i = 1:N
    for i = 1:1
         
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        [sample,error_flag] = par_aod_retri(Date,Path,Orbit,Block,r,Method,core,const,0.05);
        
        if (error_flag == 1)
            fprintf(strcat(Date,'_P',num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_B',num2str(Block),' is skipped!\n'))
            continue
        else
            save2cache(Date,Path,Orbit,Block,r,sample,Method)
        end
        
    end

end