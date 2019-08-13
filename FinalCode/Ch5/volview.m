function varargout = volview(I,isovalue,clr)

%[isovalue,varargin]  = getopt('isovalue',0,varargin{:});
%[fig,varargin]       = getopt('fig',gcf,varargin{:});
%[facecolor,varargin] = getopt('facecolor',.75*[1,1,1],varargin{:});
%[facealpha,varargin] = getopt('facealpha',.99,varargin{:});

%isovalue 
fig = 1;


I(I<isovalue) = 0;

figure(fig);
p(1) = patch(isosurface(I,isovalue));

% switch computer,
% case 'PCWIN1',
%   set(p(1),...
%     'FaceColor',facecolor,...
%     'EdgeColor','none');
% otherwise,
%   set(p(1),...
%     'FaceColor',facecolor,...
%     'EdgeColor','none',...
%     'FaceAlpha',facealpha);
% end;

isonormals(I,p(1))

%p(2) = patch(isocaps(I,isovalue),'FaceColor', 'red', 'EdgeColor', 'none');
p(2) = patch(isosurface(I,isovalue),'FaceColor', clr, 'EdgeColor', clr);

%set(p(2),...
%  'FaceColor','interp',...
%  'EdgeColor','none');

view(3); 
axis equal
axis([1,size(I,2),1,size(I,1),1,size(I,3)]);
camlight; lighting phong
material dull
drawnow

if nargout>0,
  varargout{1} = p;
end;


function [erg,arg] = getopt(str,default,varargin)
%function [erg,arg] = getopt(str,default,varargin)
% JM: 2004/02/25
% [var,varargin] = getopt(var,default,varargin{:});
% sets erg -> default, 
% if str = varargin{j}, erg = varargin{j+1}, varargin{j,j+1} = {};

erg = default;
arg = varargin;

j   = max(find(strcmp(varargin,str)));
if ~isempty(j),
  erg = varargin{j+1};
  if j == length(varargin),
    arg = {varargin{1:j-1}};
  else
    arg = {varargin{1:j-1},varargin{j+2:end}};
  end;
end;

%-------------------------------------------------------------------

