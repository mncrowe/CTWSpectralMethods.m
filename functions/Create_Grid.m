function grid = Create_Grid(H,Nz,type_z,Ly,Ny,type_y,Hy)
% Creates a grid structure containing of the coordinate grids y and z. The
% structure also contains the coordinate transformation variables zeta and 
% lambda as well as the spectral differentiation matrices for the 
% coordinate transformed grid. The offshore (y) grid may consist of an 
% arbitrary number of segments.
%
% Inputs:
%
% - H: function or constant, H = H(y) or H = H0, describing the shelf depth
% - Nz: number of gridpoints in vertical domain
% - type_z: type of vertical grid, recommended Chebyshev of second kind = 2
% - Ly,Ny: horizontal widths and number of gridpoints in each segment
% - type_y: type of horizontal grid segments, recommended Chebyshev of
%           second kind (2) or Laguerre (4, final segment only)
%
% -------------------------------------------------------------------------
% Notes:
%
% This function requires 'grid_spectral.m' and 'grid_composite.m' which
% create and join spectral basis and discretised differentiation operators.
% If you want to create custom bases, use these functions directly. A
% composite grid (i.e. one consisting of several segments) must have
% endpoints that are consistent, i.e. joining y1 and y2 required y1(end) =
% y2(1), this is only true for certain bases as some do not inculde the
% endpoints. As such, only type 2 and 4 will work in general.

if isfloat(H) == 1; H = @(y) 0*y+H; end

[Mzeta, zeta] = grid_spectral(type_z,Nz,[-1 0]);

if length(Ly) == length(type_y)
    Ly = [0 Ly];
end

for i = 1:length(Ly)-1

    if i == 1 && type_y(i) == 4 && length(type_y) > 1
        flip = 1; % flip Laguerre segment if they're the first segment of a multi-segment grid
        L = [sum(Ly(1:i+1)) sum(Ly(1:i))];
    else
        flip = 0;
        L = [sum(Ly(1:i)) sum(Ly(1:i+1))];
    end 

    [Myi, yi] = grid_spectral(type_y(i),Ny(i),L,flip);

    if i == 1
        Mlambda = Myi; lambda = yi;
    else    % join grids and differentiation matrices if there are multiple segments
        [Mlambda,lambda] = grid_composite(lambda,Mlambda,yi,Myi);
    end
end

y = lambda*ones(1,Nz);
z = H(lambda)*(zeta');

grid.y = y;
grid.z = z;
grid.zeta = zeta;
grid.lambda = lambda;
grid.Mlambda = Mlambda;
grid.Mzeta = Mzeta;
grid.H = H;
grid.Ly = Ly;
grid.type_z = type_z;
grid.type_y = type_y;
if nargin > 6; grid.Hy = Hy; else; grid.Hy = 0; end

end