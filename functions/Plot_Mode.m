function Plot_Mode(p,y,z)
% Plots the input p as a surface and contour plot in y and z.

figure;

pcolor(y,z,p)
hold on
xlabel('y'); ylabel('z')
shading interp
contour(y,z,p,[-0.8,-0.4,-0.2,0,0.2,0.4,0.8],'k','LineWidth',0.7)
colormap(cmap2([],0,[],[],0))
colorbar
line([y(1,1) y(end,1)],[z(end,1) z(end,1)],'Color','black');
line([y(1,1) y(end,1)],[z(end,end) z(end,end)],'Color','black');
line([y(1,1) y(1,1)],[z(end,1) z(end,end)],'Color','black');
line([y(end,1) y(end,1)],[z(end,1) z(end,end)],'Color','black');
hold off

set(gca,'FontSize',12,'linewidth',0.7);

end