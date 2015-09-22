function [aoda,aodb] = load_aod_batch(Location,const,Opt)
    
    [NUM,TXT,~] = xlsread('src/MISR_INFO.xls');
    
    id = find(strcmp(TXT(2:end,7),Location));

    Dates = TXT(id+1,2);
    Paths = NUM(id+1,3);
    Orbits = NUM(id+1,4);
    Blocks = NUM(id+1,5);
    
    %N = length(Dates);
    aoda = [];
    aodb = [];
    %theta=[];
    %cor1 = [];
    %cor2 = [];
    %cor3 = [];
    %ratio = [];
    %d = [];
    
    for i = [8,9,11,12,14]
        Date = Dates{i};
        Path = Paths(i);
        Orbit = Orbits(i);
        Block = Blocks(i);
        
        [aod1, ~, x1, y1, ~, ~] = load_aod(Date,Path,Orbit,Block,const,Opt);
        [aod2, x2, y2, ~, ~] = load_aeronet(Date,Path,Block,Location,const);
        
        if isempty(aod1) || isempty(aod2)
            fprintf(strcat(Date,'_P',num2str(Path,'%03d'),'_O',num2str(Orbit,'%06d'),'_B',num2str(Block),'\n',Opt,' :',num2str(length(aod1)),' aeronet:',num2str(length(aod2)),' is skipped!\n'))
            continue
        else
            [~,~,I1,I2] = match_aeronet(x1,y1,x2,y2);
            fprintf('%s:%d,aeronet:%d,%d points are matched!\n',Opt,length(x1),length(x2),length(I1))
            aoda = [aoda;aod1(I1)];
            aodb = [aodb;aod2(I2)];
            
            %theta = [theta,theta1(:,I1)];
            %[tmp1,tmp2,tmp3,tmp4] = extract_cor(reg,sample,xid,yid,const);
            %cor1 = [cor1;tmp1'];
            %cor2 = [cor2;tmp2'];
            %cor3 = [cor3;tmp3'];
            %ratio = [ratio;tmp4'];
            %tmp5 = cnt*ones(length(I1),1);
            %d = [d;tmp5];

        end
    end
    
    %varargout{1} = cor1;
    %varargout{2} = cor2;
    %varargout{3} = cor3;
    %varargout{4} = ratio;
    %varargout{5} = theta;
    %varargout{6} = d;

end
