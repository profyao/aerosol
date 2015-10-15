function plot_1d(r, x_1d, xid, yid, cmap, const, varargin)
    
    XDim_r = const.XDim_r4400 * const.r4400/r;
    YDim_r = const.YDim_r4400 * const.r4400/r;

    data = NaN * ones(XDim_r,YDim_r);
    data(sub2ind([XDim_r,YDim_r],xid,yid)) = x_1d;
    
    if nargin==7
        h = imagesc(data,varargin{1});colormap(cmap),colorbar,pbaspect([4 1 1])
    else
        h = imagesc(data);colormap(cmap),colorbar,pbaspect([4 1 1])
    end
    
    set(h, 'AlphaData', ~isnan(data))
    
end