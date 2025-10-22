function Plot_Mode(p,y,z,H,y_label,z_label,y_range,z_range)
% Plots the input p as a surface and contour plot in y and z.

arguments
    p (:,:)
    y (:,:)
    z (:,:)
    H (:,1) = 1
    y_label = 'y'
    z_label = 'z'
    y_range = [min(y,[],"all"), max(y,[],"all")];
    z_range = [min(z,[],"all"), max(z,[],"all")];
end

figure;

pcolor(y,z,p)
hold on
xlabel(y_label); ylabel(z_label)
shading interp
plot(y(:,1),-H,'k','LineWidth',0.7)
colormap(cmap2([],0,[],[],0))
cb = colorbar;
cb_tick = get(cb,"XTick");
if abs(min(p,[],"all") - cb_tick(1)) < 0.05 * abs(min(p,[],"all")); cb_tick = cb_tick(2:end); end
if abs(max(p,[],"all") - cb_tick(end)) < 0.05 * abs(max(p,[],"all")); cb_tick = cb_tick(1:end-1); end
contour(y,z,p,cb_tick,'k','LineWidth',0.7)
rectangle('Position',[y_range(1) z_range(1) y_range(2)-y_range(1) z_range(2)-z_range(1)],'LineWidth',0.7)
xlim(y_range); ylim(z_range)
hold off

set(gca,'FontSize',12,'linewidth',0.7);

end