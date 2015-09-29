function [aoda,aodb] = load_aod_rep(Location,const,Method,test_delta)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    if strcmp(Method,'MCMC')
        rep_num = 1;
    elseif test_delta==0
        rep_num = 100;
    else
        rep_num = 6;
    end
    
    aoda = [];
    aodb = [];
    
    if test_delta == 0
        load(strcat('cache/result/',Method,'_boot.mat'),'boot')
    else
        load(strcat('cache/result/',Method,'_delta.mat'),'boot')
    end
      
    cnt = 1;

    for i = [8,9,11,12,14]

        Date = Dates{i};
        Path = Paths(i);
        Block = Blocks(i);
        Orbit = Orbits(i);
        [~, ~, x1, y1, ~, ~] = load_aod(Date,Path,Orbit,Block,const,Method);
        [aod2, x2, y2, ~, ~] = load_aeronet(Date,Path,Block,Location,const);
        [~,~,I1,I2] = match_aeronet(x1,y1,x2,y2);
        fprintf('%s:%d,aeronet:%d,%d points are matched!\n',Method,length(x1),length(x2),length(I1))
        
        if strcmp(Method,'MCMC')            
            tmp = boot{cnt};
            aod1 = tmp(I1,500:5:end);
        else
            aod1 = [];
            for r = 1:rep_num
                tmp = boot{cnt,r};
                aod1 = [aod1,tmp(I1)];
            end
        end
        
        aoda = [aoda;aod1];
        aodb = [aodb;aod2(I2)];
        
        cnt = cnt + 1;
    end
end
