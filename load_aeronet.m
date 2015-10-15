function [aod, xid, yid, lon_a,lat_a] = load_aeronet(Date,Path,Block,r,Location,const) 
    
    [lon,lat] = get_coord(Path,Block,const);
    RegSize = r/const.r1100;
    file_aeronet = fullfile('aeronet/processed',Location, [Date,'_aeronet.csv']);
    
    try
        aeronet = csvread(file_aeronet,1,1);
        lon_a = aeronet(:,1);
        lat_a = aeronet(:,2);
        aod = aeronet(:,3:6);

        X = double([lon(:),lat(:)]);
        tri = delaunayn(X);
        q = [lon_a, lat_a];
        idx = dsearchn(X,tri,q);
        m = length(lon_a);
        xid = NaN*ones(m,1);
        yid = NaN*ones(m,1);

        for i=1:m
            
            if mod(idx(i),const.XDim_r1100) == 0
                xid(i) = ceil(const.XDim_r1100/RegSize);
            else
                xid(i) = ceil(mod(idx(i),const.XDim_r1100)/RegSize);
            end

            yid(i) = ceil(ceil(idx(i)/const.XDim_r1100)/RegSize);
                
        end
        
    catch
        aod = [];
        xid = [];
        yid =[];
        lon_a=[];
        lat_a=[];
    end
        

end