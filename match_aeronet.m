function [xid,yid,Ia,Ib] = match_aeronet(x0,y0,x1,y1)

    [coord,Ia,Ib] = intersect([x0,y0],[x1,y1],'rows');
    xid = coord(:,1);
    yid = coord(:,2);
end
