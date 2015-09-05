function plot_1d(x_1d, xid, yid, cmap, const, varargin)
    
    data = NaN * ones(const.XDim_r,const.YDim_r);
    data(sub2ind([const.XDim_r,const.YDim_r],xid,yid)) = x_1d;
    
    if nargin==6
        h = imagesc(data,varargin{1});colormap(cmap),colorbar,pbaspect([4 1 1])
    else
        h = imagesc(data);colormap(cmap),colorbar,pbaspect([4 1 1])
    end
    
    set(h, 'AlphaData', ~isnan(data))
    
end