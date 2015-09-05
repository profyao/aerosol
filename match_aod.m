function [aoda,aodb,xid,yid,Ia,Ib] = match_aod(aod0,x0,y0,aod1,x1,y1)

    [coord,Ia,Ib] = intersect([x0,y0],[x1,y1],'rows');
    aoda = aod0(Ia);
    aodb = aod1(Ib);
    xid = coord(:,1);
    yid = coord(:,2);
end
