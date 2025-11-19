function h = color_line3(x, y, z, c, varargin)
% color_line3 plots a 3-D "line" with c-data as color
%
%       h = color_line3(x, y, z, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line3(x, y, z, c, mark) 
%          or
%       h = color_line3(x, y, z, c, 'Colormap', cmap, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x        x-data
%       y        y-data
%       z        z-data
%       c        4th dimension for colouring
%       mark     for scatter plots with no connecting line
%       Colormap optional custom colormap (Nx3 array)
%
% out:  h        handle of the surface object

% Default colormap (jet)
cmap = parula(256);

% Check for 'Colormap' argument
idx = find(strcmpi(varargin, 'colormap'), 1);
if ~isempty(idx)
    cmap = varargin{idx + 1};
    varargin([idx idx+1]) = []; % Remove 'Colormap' and its value
end

% Normalize c values to [1, size(cmap,1)]
c = c(:);
c_scaled = round(rescale(c, 1, size(cmap,1)));
c_colors = cmap(c_scaled, :);

% Create surface plot with colors
h = surface(...
    'XData',[x(:) x(:)],...
    'YData',[y(:) y(:)],...
    'ZData',[z(:) z(:)],...
    'CData',[c(:) c(:)],...
    'FaceColor','none',...
    'EdgeColor','flat',...
    'Marker','none');

% Apply colormap to the current figure
colormap(cmap);

% Handle marker or other properties
if ~isempty(varargin)
    if ischar(varargin{1}) && ismember(varargin{1}, {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'})
        set(h, 'LineStyle', 'none', 'Marker', varargin{1})
    else
        set(h, varargin{:})
    end
end
