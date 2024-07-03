% Quasi-barotropic instabilities of a coastal jet

close all
clear
addpath("functions","Barotropic_Instability")

% define grid and parameters:
Ny = [21 31];
Nz = 21;
Ly = [1 7];
H = @(y,H0) H0 + (1-H0)*tanh(y);
U = @(y,z) erf(2*(y-1));
n = 6;

% plots dispersion curves for range of H_0 values:
plot_disp = 1;

if plot_disp == 1

    grid_baro = Create_Grid_Barotropic(Ly,Ny,[2 4]);
    params_baro = Create_Params_Barotropic(@(y) U(y,0),0,0,0);
    
    % find barotropic speed/growth rate curves as function of k:
    k = 0:0.05:1.5;
    omega0 = 0*k;
    
    for ik = 2:length(k)
        omega_t = Find_Modes_Barotropic(grid_baro,params_baro,k(ik),n,1i/pi,0,'lm');
        omega0(ik) = omega_t(imag(omega_t)==max(imag(omega_t)));
    end
    
    % plot barotropic curves:
    figure; plot(k,real(omega0),k,imag(omega0),'LineWidth',0.7); grid on
    xlabel('k'); legend('Re[\omega]','Im[\omega]')
    set(gca,'FontSize',12,'linewidth',0.7);
    
    % find closest 2D modal solutions:
    H0 = 1:-0.1:0.1;
    omega = [omega0.' zeros(length(k),length(H0)-1)];
    
    params_yz = Create_Params(U,@(z) 1+0*z,1,1,0,0);
    for iH = 2:length(H0)
    
        grid_yz = Create_Grid(@(y) H(y,H0(iH)),Nz,2,Ly,Ny,[2 4]);
    
        for ik = 2:length(k)
    
            disp([' - Finding solution for (k,H_0) = (' num2str(k(ik)) ',' num2str(H0(iH)) ')'])
    
            % determine initial guess for omega using interpolation:
            if iH > 2 && ik > 2
                w0 = 0.5*(2*omega(ik,iH-1) - omega(ik,iH-2)) + 0.5*(2*omega(ik-1,iH) - omega(ik-2,iH));
            else
                w0 = omega(ik,iH-1);
            end
    
            [omega_t,p,~,~,~,~,~,~,~,~,~,~,zero_wall] = Find_Modes(grid_yz,params_yz,k(ik),n,w0,0,'lm');
            [~,i] = max(sum(abs(p(1,:,:)/H0(iH))));    % quasi-barotropic mode has pressure signal at y = 0                                                     
            %[~,i] = min(abs(omega_t - omega(ik,iH-1)));    % alternative method, find closest omega to guess
            omega(ik,iH) = omega_t(i);
    
        end
    end
    
    % plot solution curves for H_0 = 1, 0.7, 0.4, 0.1:
    figure
    plot(k,real(omega(:,[1 4 7 10])),'LineWidth',0.7)
    xlabel('k'); ylabel('Re[\omega]'); grid
    legend('H_0 = 1.0','H_0 = 0.7','H_0 = 0.4','H_0 = 0.1')
    set(gca,'FontSize',12,'linewidth',0.7);
    
    figure
    plot(k,imag(omega(:,[1 4 7 10])),'LineWidth',0.7)
    xlabel('k'); ylabel('Im[\omega]'); grid
    legend('H_0 = 1.0','H_0 = 0.7','H_0 = 0.4','H_0 = 0.1')
    set(gca,'FontSize',12,'linewidth',0.7);

end

% plot modes for one example setup:
plot_mode = 1;          % 1 - plot specified mode
H0_plot = 0.5;          % value of H0 to plot
k_plot = 1;             % wavenumber to plot
Nz_plot = 31;
Ny_plot = [31 51];

[~,iH] = min(abs(H0-H0_plot));
[~,ik] = min(abs(k-k_plot));
w0 = omega(ik,iH);      % initial guess for omega, use previous results

if plot_mode == 1
    
    grid_yz = Create_Grid(@(y) H(y,H0_plot),Nz_plot,2,Ly,Ny_plot,[2 4]);
    params_yz = Create_Params(U,@(z) 1+0*z,1,1,0,0);

    [omega_plot,p,u,v,w,b,U0,U0_y,U0_z,~,y,z] = Find_Modes(grid_yz,params_yz,k_plot,n,w0,0,'lm');
    [~,plot_mode] = min(abs(omega_plot-w0));
    
    Plot_Mode(abs(p(:,:,plot_mode))/max(max(abs(p(:,:,plot_mode)))),y,z)
    disp(['Plot of |p| for (k,H_0) = (' num2str(k_plot) ',' num2str(H0_plot) ')'])
end