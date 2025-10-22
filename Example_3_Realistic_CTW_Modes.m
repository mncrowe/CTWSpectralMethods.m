% Display realistic CTW modes based on Gelderloos et. al. (2021) setup

close all
clear
addpath("functions")

% set parameters:
f = 1e-4;
g = 10;
Nz = 21;
Ly = [1e5 2e5];
Ny = [31 31];
type_y = [2 4];

% unit conversion (time):
cyc_day_to_rad_s = @(w) w*(2*pi)/60/60/24;
rad_s_to_cyc_day = @(w) w/(2*pi)*60*60*24;

% unit conversion (space):
cyc_km_to_rad_m = @(k) k*(2*pi)/1e3;
rad_m_to_cyc_km = @(k) k/(2*pi)*1e3;

% N2 from exponential fit to Gelderloos data:
N2 = @(y,z) 2e-4*exp(z/90)+1e-6;

% H from smoothed, windowed profile from Gelderloos data:
load('Gelderloos_et_al_dat.mat')
y1 = 0.08e5; y2 = 1.85e5;   % use H profile for y in [y1,y2], otherwise flat bottom
H = @(y) smooth((ones(size(y))*interp1(yd,Hd,y1)).*(y<=y1)+interp1(yd,Hd,y).*(y<y2).*(y>y1)+(ones(size(y))*interp1(yd,Hd,y2)).*(y>=y2));

% U from exponential/Gaussian fit to Gelderloos data:
U = @(y,z) 0.4*exp(-(y-9.5e4).^2/3e4^2).*exp(z/200);

% M2 from U via thermal wind (M2 = -f*U_z):
M2 = @(y,z) -f * 0.4*exp(-(y-9.5e4).^2/3e4^2).*exp(z/200) / 200;

% create grid and parameters, use H between y = 1e4 and y = 1.85e5:
grid = Create_Grid(H,Nz,2,Ly,Ny,type_y);
params = Create_Params(U,M2,N2,f,g,1,1);

% plot stratification:
figure
z_plot = -(0:10:1800);
semilogx(N2(0,z_plot),z_plot,'k','LineWidth',0.7)
grid on
xlabel('N^2 (s^{-2})'); ylabel('z (m)')
xlim([7e-7 2e-4])
set(gca,'FontSize',12,'linewidth',0.7);

% plot velocity:
Plot_Mode(U(grid.y,grid.z),grid.y/1e3,grid.z,H(grid.y(:,1)),'y (km)','z (m)',[0 2e2])
cb = colorbar; ylabel(cb,'U (m/s)')

% create wavenumber domain and initial guesses for each mode:
k = cyc_km_to_rad_m((0.2:0.2:7)*1e-3);
omega = zeros(length(k),3);
p = zeros(length(k),3,sum(Ny)-length(Ny)+1,Nz);
M = 3;      % number of modes to find
n = 6;

% initial values (determined by examining all solutions for k = 0.2e-3):
omega(1,1) = cyc_day_to_rad_s(0.093903153350806);
omega(1,2) = cyc_day_to_rad_s(0.015539651145854);
omega(1,3) = cyc_day_to_rad_s(0.008657303700164);

% loop through modes and k values to find frequency curves for each mode:
for im = 1:M
    
    disp(['- Mode ' num2str(im) ', finding omega for k = ' num2str(rad_m_to_cyc_km(k(1)))])
    [omega_t,p_t] = Find_Modes(grid,params,k(1),n,omega(1,im),0,'lm');
    [~,i] = min(abs(omega_t - omega(1,im)));
    omega(1,im) = omega_t(i);
    p(1,im,:,:) = p_t(:,:,i);

    for ik = 2:length(k)

        disp(['- Mode ' num2str(im) ', finding omega for k = ' num2str(rad_m_to_cyc_km(k(ik)))])
        
        % estimate frequency to use as initial guess:
        if ik == 2
            w0 = omega(ik-1,im)*k(ik)/k(ik-1);
        else
            w0 = 2*omega(ik-1,im)-omega(ik-2,im);
        end

        [omega_t,p_t] = Find_Modes(grid,params,k(ik),n,w0,0,'lm');

        % remove (numerical) instabilities:
        i_instab = (imag(omega_t) == 0);
        omega_t = omega_t(i_instab); p_t = p_t(:,:,i_instab);

        % limit difference in pressure eigenfunctions:
        p_diff = squeeze(sum(abs(p_t - squeeze(p(ik-1,im,:,:))),[1 2]))/(Nz*(sum(Ny)-length(Ny)+1));
        i_pdiff = (p_diff < 0.15);
        omega_t = omega_t(i_pdiff); p_t = p_t(:,:,i_pdiff);

        % find closest frequency to guess and save:
        [~,i] = min(abs(omega_t-w0));
        omega(ik,im) = omega_t(i);
        p(ik,im,:,:) = p_t(:,:,i);

    end

end

% plot omega against k for calculated modes:
figure
plot(rad_m_to_cyc_km(k),rad_s_to_cyc_day(omega),'LineWidth',0.7); grid on
xlabel('k (cyc/km)'); ylabel('\omega (cyc/day)')
legend('Mode I','Mode II', 'Mode III')
set(gca,'FontSize',12,'linewidth',0.7);

% find a single Mode III solution:
k_plot = cyc_km_to_rad_m(5*1e-3);
[~,p_plot] = Find_Modes(grid,params,k_plot,n,omega(k==k_plot,3),0,'lm');

% plot the pressure eigenfunction for this Mode III solution:
Plot_Mode(p_plot(:,:,1),grid.y/1e3,grid.z,H(grid.y(:,1)),'y (km)','z (m)',[0 2e2])