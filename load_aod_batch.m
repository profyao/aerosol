function [aoda,aodb] = load_aod_batch(Method,Location,const)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    N = length(Dates);
    aoda = [];
    aodb = [];
    
    for i = 1:N
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        [aod1, x1, y1, ~, ~] = load_aod(Date,Path,Orbit,Block,const,Method);
        [aod2, x2, y2, ~, ~] = load_aeronet(Date,Path,Block,Location,const);
        
        if isempty(aod1) || isempty(aod2)
            fprintf(strcat(Date,'_P',num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_B',num2str(Block),'\n',Method,' :',num2str(length(aod1)),' aeronet:',num2str(length(aod2)),' is skipped!\n'))
            continue
        else
            [tmpa,tmpb,~,~,~,~] = match_aod(aod1,x1,y1,aod2,x2,y2);
            fprintf('%s:%d,aeronet:%d,%d points are matched!\n',Method,length(aod1),length(aod2),length(tmpa))
            aoda = [aoda;tmpa];
            aodb = [aodb;tmpb];
        end
    end

end
