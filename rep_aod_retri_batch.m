function rep_aod_retri_batch(Method,r,Location,core,test_delta,const)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));   

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    if strcmp(Method,'MCMC')
        rep_num = 1;
    elseif test_delta==0
        rep_num = 1000;
    else
        delta_all = [0.001,0.01,0.05,0.1,1,10];
        rep_num = 6;
    end

    boot = cell(5,rep_num);
    
    for rep = 1:rep_num
        
        cnt = 1;
        if test_delta == false
            delta = 0.05;
        else
            delta = delta_all(rep);
            fprintf('delta is %f \n',delta)
        end

        for i = [8,9,11,12,14]
            
            Date = Dates{i};
            Path = Paths(i);
            Orbit = Orbits(i);
            Block = Blocks(i);

            disp([cnt,rep])
            [sample,~] = par_aod_retri(Date,Path,Orbit,Block,r,Method,core,const,delta);
            boot{cnt,rep} = sample.tau;
            clc
            
            cnt = cnt + 1;
        end
    
    end
    
    if test_delta == 0
        save(strcat('cache/result/',Method,'_boot.mat'),'boot')
    else
        save(strcat('cache/result/',Method,'_delta.mat'),'boot')
    end

end