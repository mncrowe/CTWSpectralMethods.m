function [omega,p,u,v,U,Uy,y] = Find_Modes_Barotropic(grid,params,k,n,w0,BC,method)
% Solves the barotropic wave/instability problem
%
% Inputs:
%
% - grid: structure created using 'Create_Grid_Barotropic.m' containing grid variables
% - params: structure created using 'Create_Params_Barotropic.m' containing physical parameters
% - k: along-shore (x) wavenumber
% - n: number of modes to find
% - w0: initial guess for frequency
% - BC: boundary conditions at y = y_max, 0: none (default), 1: v = 0, 2: v_y = 0
% - method: string, method for Gen_EVP (default: 'lm')
%
% Outputs:
%
% - omega: vector of frequencies (inc. growth rates if complex)
% - (p,u,v): pressure, velocity and buoyancy fields for each mode
% - (U,Uy): background flow plus gradient
% - y: offshore coordinate vector

if nargin < 3; k = 1; end
if nargin < 4; n = 10; end
if nargin < 5; w0 = 1/pi; end
if nargin < 6; BC = 0; end
if nargin < 7; method = 'lm'; end

y = grid.y; My = grid.My;

N = length(y);
U = params.U(y);
Uy = My*(U-U(end));
f = params.f;
nu = params.nu;

I = eye(N);
O = zeros(N);

L = [I O O; O -I O; O O O];
R = [k*diag(U)+1i*nu*My^2 diag(Uy-f) k*I; f*I -k*diag(U)-1i*nu*My^2 My; k*I My O];

L(N+1,:) = 0;
R(N+1,:) = [zeros(1,N) 1 zeros(1,2*N-1)];

if nu > 0
    L(1, :) = [My(1, :) zeros(1, 2*N)]; % no stress on u
    R(1, :) = 0;
end

if BC == 1
    R(end,:) = [zeros(1,2*N-1) 1 zeros(1,N)];
    L(end,:) = 0;
end

if BC == 2
    R(end,:) = [zeros(1,N) My(end,:) zeros(1,N)];
    L(end,:) = 0;
end

[omega,phi] = Gen_EVP(L,-R,n,w0,method);

phi = reshape(phi,N,3,[]);
phi = phi./max(abs(phi(:,3,:))); % normalise fields

u = squeeze(phi(:,1,:));
v = 1i*squeeze(phi(:,2,:));   % system is solved for i*v instead of v so matrices are real
p = squeeze(phi(:,3,:));

end