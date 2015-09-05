function plot_2d(x_2d,cmap,varargin)
        
    if nargin==3
        h = imagesc(x_2d,varargin{1});colormap(cmap),colorbar,pbaspect([4 1 1])
    else
        h = imagesc(x_2d);colormap(cmap),colorbar,pbaspect([4 1 1])
    end
    
    set(h, 'AlphaData', ~isnan(x_2d))
end