function [lon,lat] = conv_coord(lon,lat,res_down_multiple,const)

    c = res_down_multiple;
    
    if c > 1
        
        left = round(c/2);
        right = left+1;
    
        lon1 = lon(left:c:const.XDim_r1100,left:c:const.YDim_r1100);
        lon2 = lon(left:c:const.XDim_r1100,right:c:const.YDim_r1100);
        lon3 = lon(right:c:const.XDim_r1100,left:c:const.YDim_r1100);
        lon4 = lon(right:c:const.XDim_r1100,right:c:const.YDim_r1100);
        lon = (lon1+lon2+lon3+lon4)*0.25;

        lat1 = lat(left:c:const.XDim_r1100,left:c:const.YDim_r1100);
        lat2 = lat(left:c:const.XDim_r1100,right:c:const.YDim_r1100);
        lat3 = lat(right:c:const.XDim_r1100,left:c:const.YDim_r1100);
        lat4 = lat(right:c:const.XDim_r1100,right:c:const.YDim_r1100);
        lat = (lat1+lat2+lat3+lat4)*0.25;
        
    elseif c < 1
        
        error('resolution not implemented!')
        
    end

end