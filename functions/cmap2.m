function map = cmap2(R,centre,f,N,rescale,C,debug)
% outputs a custom colour map for use with colormap(map)
%
% - R: plot limits, integer or array;
%    integer: use current figure CLim value (default)
%    vector: 2 element vector of [Min Max]
%    field: F, calculate limits as min and max of F array
% - centre: value correspnding to the colorbar centre (default: (Min+Max)/2)
% - f: function f:[0,1]->[0,1], nonlinear rescaling (default: f = @(x) x)
% - N: size of colorbar (default: 256)
% - rescale: 1 - rescale ends to max/min colour bar values, 0 - no scaling (default: 1)
% - C: colours, integer or array;
%    integer: 0 - blue/white/red (default)
%             1 - greyscale
%    array: list of N colours (rgb), N x 3, values in [0,1], can use Matlab presets (e.g. jet(8), parula(3))
% - debug: 1 - outputs limits and centre value, 0 - no output (default)
%
% Note: default parameters may be entered through undefined variables or by
% setting the variable to [] in input. Requires dist_array.m and F_minmax.m

% define unspecified variables
if nargin < 1; R = []; end
if nargin < 2; centre = []; end
if nargin < 3; f = []; end
if nargin < 4; N = []; end
if nargin < 5; rescale = []; end
if nargin < 6; C = []; end
if nargin < 7; debug = 0; end

% set variables defined using flags and input arrays
if numel(R) < 2; [Min,Max] = dist_array(get(gca,'CLim')); end
if numel(R) == 2; Min = R(1); Max = R(2); end
if numel(R) > 2; [Min,Max] = F_minmax(R); end
if numel(centre) == 0; centre = (Min+Max)/2; end
if isa(f,'function_handle') == 0; f = @(x) x; end
if numel(N) == 0 || N == 0; N = 256; end
if numel(rescale) == 0; rescale = 1; end
if numel(C) == 0; C = 0; end
if numel(C) == 1
    if C == 0; C = [0 0 0.5; 0 0.5 1; 1 1 1; 1 0 0; 0.5 0 0];  end
    if C == 1; C = [1 1 1; 0.5 0.5 0.5; 0 0 0]; end
end

Tick_lims = (linspace(Min,Max,N)-centre)/max(Max-centre,centre-Min);
if rescale == 1    % rescale ends if rescale == 1
    if Tick_lims(1) > -1; Tick_lims(Tick_lims<0) = -Tick_lims(Tick_lims<0)/Tick_lims(1); end
    if Tick_lims(end) < 1; Tick_lims(Tick_lims>0) = Tick_lims(Tick_lims>0)/Tick_lims(end); end
end

map = interp1(linspace(-1,1,length(C)),C,sign(Tick_lims).*f(abs(Tick_lims)));   % create colour map

if debug == 1   % debugging, print outputs to screen
    disp(['Limits = [' num2str([Min Max]) ']'])
    disp(['Centre = ' num2str(centre)])
end

end

