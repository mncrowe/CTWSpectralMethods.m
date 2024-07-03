function [omega,p,u,v,w,b,U,Uy,y] = Find_Modes_Equiv_Barotropic(grid,params,k,m,n,w0,BC,method)
% Solves the equivalent barotropic wave/instability problem
%
% Inputs:
%
% - grid: structure created using 'Create_Grid.m' containing grid variables
% - params: structure created using 'Create_Params.m' containing physical parameters
% - k: along-shore (x) wavenumber
% - m: vertical wavenumber, m = M/H for integer M
% - n: number of modes to find
% - w0: initial guess for frequency
% - BC: boundary conditions at y = y_max, 0: none (default), 1: v = 0, 2: v_y = 0
% - method: string, method for Gen_EVP (default: 'lm')
%
% Outputs:
%
% - omega: vector of frequencies (inc. growth rates if complex)
% - (p,u,v,w,b): pressure, velocity and buoyancy fields for each mode
% - (U,Uy): background flow plus gradient
% - y: offshore coordinate vector
%
% Note: the full z dependent fields can be recovered using:
%
%   u = u'(y) cos(m pi z) exp(i(k x - omega t))
%   v = v'(y) cos(m pi z) exp(i(k x - omega t))
%   w = w'(y) sin(m pi z) exp(i(k x - omega t))
%   b = b'(y) sin(m pi z) exp(i(k x - omega t))
%   p = p'(y) cos(m pi z) exp(i(k x - omega t))
%
% where (u',v',w',b',p') are the outputs from this function. The value of m
% should be chosen such that sin(m pi H) = 0. i.e. m = M/H where M is an
% integer.

if nargin < 3; k = 1; end
if nargin < 4; m = 1; end
if nargin < 5; n = 10; end
if nargin < 6; w0 = 1/pi; end
if nargin < 7; BC = 0; end
if nargin < 8; method = 'lm'; end

y = grid.y; My = grid.My;

N = length(y);
U = params.U(y);
Uy = My*(U-U(end));
f = params.f;
N2 = params.N2;
if params.hydrostatic; h = 1; else; h = 0; end

I = eye(N);
O = zeros(N);

L = [I O O O O; O -I O O O; O O -h*I O O; O O O I O; O O O O O];
R = [k*diag(U) diag(Uy-f) O O k*I; f*I -k*diag(U) O O My; O O -h*k*diag(U) -I -m*pi*I; O O N2*I k*diag(U) O; k*I My m*pi*I O O];

L(N+1,:) = 0;
R(N+1,:) = [zeros(1,N) 1 zeros(1,4*N-1)];

if BC == 1
    R(end,:) = [zeros(1,2*N-1) 1 zeros(3,N)];
end

if BC == 2
    R(end,:) = [zeros(1,N) My(end,:) zeros(3,N)];
end

[omega,phi] = Gen_EVP(L,-R,n,w0,method);

phi = reshape(phi,N,5,[]);
phi = phi./max(abs(phi(:,5,:))); % normalise fields

u = squeeze(phi(:,1,:));
v = 1i*squeeze(phi(:,2,:));   % system is solved for i*v and i*w instead of v and w so matrices are real
w = 1i*squeeze(phi(:,3,:));
b = squeeze(phi(:,4,:));
p = squeeze(phi(:,5,:));

end