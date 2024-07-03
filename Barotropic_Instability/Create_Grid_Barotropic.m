function grid = Create_Grid_Barotropic(L,N,type)
% Creates a 1d y grid for the barotropic instability problem

if nargin < 3; type = [2 4]; end
if nargin < 2; N = [21 31]; end
if nargin < 1; L = [2 3]; end

L = [0 L];
for i = 1:length(L)-1
    [Mi, yi] = grid_spectral(type(i),N(i),[sum(L(1:i)) sum(L(1:i+1))]);
    if i == 1
        My = Mi; y = yi;
    else
        [My,y] = grid_composite(y,My,yi,Mi);
    end
end

grid.y = y;
grid.My = My;
grid.Ly = L;
grid.type = type;

end