function [lon,lat] = conv_coord(lon,lat,res_down_multiple,const)

    c = res_down_multiple;
    
    lon1 = lon(2:c:const.XDim_r1100,2:c:const.YDim_r1100);
    lon2 = lon(2:c:const.XDim_r1100,3:c:const.YDim_r1100);
    lon3 = lon(3:c:const.XDim_r1100,2:c:const.YDim_r1100);
    lon4 = lon(3:c:const.XDim_r1100,3:c:const.YDim_r1100);
    lon = (lon1+lon2+lon3+lon4)*0.25;
    
    lat1 = lat(2:c:const.XDim_r1100,2:c:const.YDim_r1100);
    lat2 = lat(2:c:const.XDim_r1100,3:c:const.YDim_r1100);
    lat3 = lat(3:c:const.XDim_r1100,2:c:const.YDim_r1100);
    lat4 = lat(3:c:const.XDim_r1100,3:c:const.YDim_r1100);
    lat = (lat1+lat2+lat3+lat4)*0.25;

end