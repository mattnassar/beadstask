function h = surf2d(data,x,y)
% function h = surf2d(data,x,y)
% 
% surf2d(data) should produce output that looks like imagesc(data)
% advantage over imagesc: saved figures are not 'rasterized'
%
% inputs:
%    x ... (optional) must have same length as cols in data
%    y ... (optional) must have same length as rows in data
%
% output:
%    h ... handle to surf object
%
% 2012/09/11 yl

% MN changed the default 'view' to [0 90] on 5-6-13.
% also added a brighter colormap, courtesy of john burke..
colormap(myJet)


% get the dimensions for surf
% note: surf defines the vertices, whereas our data defines the surfaces
%       in between -- hence the +1
nx = size(data,2)+1;
ny = size(data,1)+1;

% must do some fiddly things to deal with edge effects of surf
if ~exist('x','var') || isempty(x)
    x = -.5+1:nx;
else
    if length(x) ~= nx-1
        error('x dim does not match');
    end
    
    dx = unique(diff(x));
    if length(dx) > 1
        warning('x scale not linear. Image may look funny.');
        dx = min(dx);
    end
    x = -.5*dx + [x(:); x(end)+dx];
end

if ~exist('y','var') || isempty(y)
    y = -.5 + 1:ny;
else
    if length(y) ~= ny-1
        error('y dim does not match');
    end
    
    dy = unique(diff(y));
    if length(dy) > 1
        warning('y scale not linear. Image may look funny.');
        dy = min(dy);
    end
    y = -.5*dy + [y(:); y(end)+dy];
end

% ok, now actually plot these
h = surf(x,y,zeros(ny,nx),data,'edgecolor','none');
view([0 -90]);
set(gca,'box','on','layer','top','xgrid','off','ygrid','off', 'view', [0 90]);
axis tight;
