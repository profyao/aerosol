function imagesc2 ( img_data, varargin)
% a wrapper for imagesc, with some formatting going on for nans

% plotting data. Removing and scaling axes (this is for image plotting)
if nargin==3
    imagesc(img_data,varargin{1});colormap(varargin{2})
else
    imagesc(img_data);colormap(varargin{2})
end
