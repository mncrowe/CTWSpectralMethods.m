function grid = Create_Grid(H,Nz,type_z,Ly,Ny,type_y)
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

Ly = [0 Ly];
for i = 1:length(Ly)-1
    [Myi, yi] = grid_spectral(type_y(i),Ny(i),[sum(Ly(1:i)) sum(Ly(1:i+1))]);
    if i == 1
        Mlambda = Myi; lambda = yi;
    else
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

end