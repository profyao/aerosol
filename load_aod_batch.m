function [aoda,aodb] = load_aod_batch(Location,r,const,Opt)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    %N = length(Dates);
    aoda = [];
    aodb = [];

    
    for i = [8,9,11,12,14]
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        [aod1, ~, x1, y1, ~, ~] = load_aod(Date,Path,Orbit,Block,r,const,Opt);
        [aod2, x2, y2, ~, ~] = load_aeronet(Date,Path,Block,r,Location,const);
        
        if isempty(aod1) || isempty(aod2)
            fprintf(strcat(Date,'_P',num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_B',num2str(Block),'\n',Opt,' :',num2str(length(aod1)),' aeronet:',num2str(length(aod2)),' is skipped!\n'))
            continue
        else
            [~,~,I1,I2] = match_aeronet(x1,y1,x2,y2);
            fprintf('%s:%d,aeronet:%d,%d points are matched!\n',Opt,length(x1),length(x2),length(I1))
            aoda = [aoda;aod1(I1,:)];
            aodb = [aodb;aod2(I2,:)];
        end
    end

end
