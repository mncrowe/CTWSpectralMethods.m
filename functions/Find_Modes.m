function [omega,p,u,v,w,b,U,Uy,Uz,Nsqr,Msqr,y,z,zero_wall,zero_surf] = Find_Modes(grid,params,k,n,w0,BC,method)
% Solve for the Coastal wave/instability modes for the given grid and
% physical parameters
%
% Inputs:
%
% - grid: structure created using 'Create_Grid.m' containing grid variables
% - params: structure created using 'Create_Params.m' containing physical parameters
% - k: along-shore (x) wavenumber
% - n: number of modes to find
% - w0: initial guess for frequency
% - BC: boundary conditions at y = y_max, 0: none (default), 1: v = 0, 2: v_y = 0
% - method: string, method for Gen_EVP (default: 'lm')
%
% Outputs:
%
% - omega: vector of frequencies (inc. growth rates if complex)
% - (p,u,v,w,b): pressure, velocity and buoyancy fields for each mode
% - (U,Uy,Uz): background flow plus gradients on (y,z) grid
% - Nsqr: value of N^2 in (y,z) grid
% - (y,z): coordinate grid
% - zero_wall: number of zero intersections along the side and bottom boundary
% - zero_surf: number of zero intersections on top surface

f = params.f;
nu = params.nu;

if nargin < 3; k = 1; end
if nargin < 4; n = 10; end
if nargin < 5; w0 = f/pi; end
if nargin < 6; BC = 0; end
if nargin < 7; method = 'lm'; end

if length(BC) == 1; BC = [1 BC]; end % wall BC on left by default

% calculate derived parameters:

Nz = length(grid.zeta); Ny = length(grid.lambda);
y = grid.y; z = grid.z;

H = grid.H(grid.lambda)+1e-8*(tanh(3*grid.lambda/grid.lambda(end))-1); % plus small regularisation term

if ~isa(grid.Hy,"function_handle") % if Hy is given in 'grid', calculate it directly, otherwise evaluate by differentiating H
    Hy = grid.Mlambda*(H-H(end));
else
    Hy = grid.Hy(grid.lambda);
end

U = params.U(y,z);

if ~isa(params.Uy,"function_handle") % if Uy is given in 'params', calculate it directly, otherwise evaluate by differentiating U
    U0 = params.U(ones(Ny,1)*y(end,:),z); % U0 = U(y=L,z), subtract from U for calculating Uy on Laguerre grid
    Uy = (grid.Mlambda*(U-U0))-(grid.zeta.*(grid.Mzeta*(U-U0)').*(Hy./H)')';
else
    Uy = params.Uy(y,z);
end

if ~isa(params.Uz,"function_handle") % if Uz is given in 'params', calculate it directly, otherwise evaluate by differentiating U
    Uz = ((grid.Mzeta*U')./H')';
else
    Uz = params.Uz(y,z);
end

Nsqr = params.N2(y,z);
Msqr = params.M2(y,z);
if params.hydrostatic; h = 1; else; h = 0; end
if params.free_surface; a = 1; else; a = 0; end

% define matrix terms:

dy = create_operator(grid.Mlambda,0,eye(Nz),0) - create_operator(diag(Hy./H),0,diag(grid.zeta)*grid.Mzeta,0);
dz = create_operator(diag(1./H),0,grid.Mzeta,0);
O = zeros(Ny*Nz);
I = eye(Ny*Nz);
U1 = diag(reshape(U,[Ny*Nz 1]));
U2 = diag(reshape(Uy,[Ny*Nz 1]));
U3 = diag(reshape(Uz,[Ny*Nz 1]));
N2 = diag(reshape(Nsqr,[Ny*Nz 1]));
M2 = diag(reshape(Msqr,[Ny*Nz 1]));

% define diffusive term:

dlambda = kron(eye(Nz),grid.Mlambda);
diff = nu*dlambda^2; %nu*dy^2;

% build matrices for system omega*L*psi = R*psi:

L = kron(diag([1 -1 -h 1 0]),I);
R = [k*U1+1i*diff U2-f*I U3 O k*I; f*I -k*U1-1i*diff O O dy; O O -h*(k*U1+1i*diff) -I dz; O M2 N2 k*U1+1i*diff O; k*I dy dz O O];

% apply BCs:

% no stress/flux BCs on y = 0
if nu > 0 && BC(1) ~= 0
    iz = 1:Ny:Ny*Nz;
    dy_wall = dlambda(iz,:); %dy(iz,:); %dy(iz, :);
    
    % du/dy = 0 on y = 0
    R(iz,:) = 0;
    R(iz,1:Ny*Nz) = dy_wall;
    L(iz,:) = 0;
    
    % db/dy = 0 on y = 0
    R(3*Nz*Ny+iz,:) = 0;
    R(3*Nz*Ny+iz,(3*Nz*Ny+1):4*Ny*Nz) = dy_wall;
    L(3*Nz*Ny+iz,:) = 0;

    if h > 0
        % dw/dy = 0 on y = 0
        R(2*Nz*Ny+iz,:) = 0;
        R(2*Nz*Ny+iz,(2*Nz*Ny+1):3*Ny*Nz) = dy_wall;
        L(2*Nz*Ny+iz,:) = 0;
    end
end

% bottom BC: H_y v + w = 0 on z = -H(y)
R((2*Nz*Ny+1):(2*Nz*Ny+Ny),:) = 0;
R((2*Nz*Ny+1):(2*Nz*Ny+Ny),(Nz*Ny+1):(Nz*Ny+Ny)) = diag(Hy);
R((2*Nz*Ny+1):(2*Nz*Ny+Ny),(2*Nz*Ny+1):(2*Nz*Ny+Ny)) = eye(Ny);
L((2*Nz*Ny+1):(2*Nz*Ny+Ny),(2*Nz*Ny+1):(2*Nz*Ny+Ny)) = zeros(Ny);

% top BC: b + a*N^2/g p = 0 on z = 0
%R((5*Nz*Ny-Ny+1):(5*Nz*Ny),:) = 0;
%R((5*Nz*Ny-Ny+1):(5*Nz*Ny),(4*Nz*Ny-Ny+1):(4*Nz*Ny)) = eye(Ny);
%R((5*Nz*Ny-Ny+1):(5*Nz*Ny),(5*Nz*Ny-Ny+1):(5*Nz*Ny)) = a*params.N2(y(:,end),z(:,end))/params.g.*eye(Ny);
%L((5*Nz*Ny-Ny+1):(5*Nz*Ny),(3*Nz*Ny-Ny+1):(3*Nz*Ny)) = zeros(Ny);

% top BC: w - a/g*(-omega + kU)p = 0
R((5*Nz*Ny-Ny+1):(5*Nz*Ny),:) = 0;
R((5*Nz*Ny-Ny+1):(5*Nz*Ny),(3*Nz*Ny-Ny+1):(3*Nz*Ny)) = eye(Ny);
R((5*Nz*Ny-Ny+1):(5*Nz*Ny),(5*Nz*Ny-Ny+1):(5*Nz*Ny)) = -a/params.g.*k*U1((Nz-1)*Ny+1:Ny*Nz,(Nz-1)*Ny+1:Ny*Nz);
L((5*Nz*Ny-Ny+1):(5*Nz*Ny),(5*Nz*Ny-Ny+1):(5*Nz*Ny)) = a/params.g*eye(Ny);

if BC(1) == 1 % wall BC: v = 0  on y = 0
    R((Nz*Ny+1):Ny:(2*Nz*Ny),:) = 0;
    R((Nz*Ny+1):Ny:(2*Nz*Ny),(Nz*Ny+1):Ny:(2*Nz*Ny)) = eye(Nz);
    L((Nz*Ny+1):Ny:(2*Nz*Ny),(Nz*Ny+1):Ny:(2*Nz*Ny)) = zeros(Nz);
end

if BC(1) == 2 % wall BC: v = 0  on y = 0
    R((Nz*Ny+1):Ny:(2*Nz*Ny),:) = 0;
    R((Nz*Ny+1):Ny:(2*Nz*Ny),(Nz*Ny+1):Ny:(2*Nz*Ny)) = dy(1:Ny:Nz*Ny,:);
    L((Nz*Ny+1):Ny:(2*Nz*Ny),(Nz*Ny+1):Ny:(2*Nz*Ny)) = zeros(Nz);
end

if BC(2) == 1 % outer BC: v = 0 on y = L_max
    R((4*Nz*Ny+Ny):Ny:(5*Nz*Ny),:) = 0;
    R((4*Nz*Ny+Ny):Ny:(5*Nz*Ny),(Nz*Ny+Ny):Ny:(2*Nz*Ny)) = eye(Nz);
end

if BC(2) == 2 % outer BC: v_y = 0 on y = L_max
    R((4*Nz*Ny+Ny):Ny:(5*Nz*Ny),:) = 0;
    R((4*Nz*Ny+Ny):Ny:(5*Nz*Ny),Nz*Ny+1:2*Nz*Ny) = dy(Ny:Ny:Nz*Ny,:);
end

% find eignevalues, omega, using 'Gen_EVP.m':

[omega,phi] = Gen_EVP(L,-R,n,w0,method);

% reshape and normalise output:
if nargout > 1
    phi = reshape(phi,Ny,Nz,5,[]);
    phi = phi./max(max(abs(phi(:,:,5,:))));
    phi = phi./(phi(1,end,5,:)./abs(phi(1,end,5,:)));
end

if nargout > 1; p = squeeze(phi(:,:,5,:)); end
if nargout > 2; u = squeeze(phi(:,:,1,:)); end
if nargout > 3; v = 1i*squeeze(phi(:,:,2,:)); end  % system is solved for i*(v,w) so L is real
if nargout > 4; w = 1i*squeeze(phi(:,:,3,:)); end
if nargout > 5; b = squeeze(phi(:,:,4,:)); end

% look for number of zero intersections along side/bottom wall and surface:
if nargout > 12
    g_wall = real([squeeze(p(1,end:-1:1,:)).' squeeze(p(:,1,:)).']);
    g_surf = real(squeeze(p(:,end,:)).');
    zero_wall = sum(sign(g_wall(:,2:end).*g_wall(:,1:end-1))==-1,2);
    zero_surf = sum(sign(g_surf(:,2:end).*g_surf(:,1:end-1))==-1,2);
end

end