
%x=[0:.01:2];
%y1= exp(-x).*cos(6*pi*x);
%y2= exp(-x).*sin(6*pi*x);
%y3=-exp(-x).*cos(6*pi*x);
%y4=-exp(-x).*sin(6*pi*x);

%figure(1)
%% set figure position and size:
set(gcf,'position',[160 280 600 460])
%% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')
%% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)
%% to plot several curves on top of each other:
hold on
%% plot curves:

contour(xm,ym,psi,c);

%plot(x,y1,'-' ,'color','r'     ,'linewidth',2)
%plot(x,y2,'--','color',[0 .5 0],'linewidth',2) % [0 .5 0] is dark green in rgb
%plot(x,y3,'-.','color',[0  0 1],'linewidth',2) % [0 0 1] is light blue in rgb
%plot(x,y4,'-','color','k','linewidth',2,'markersize',8,'markerfacecolor','w')
%plot(x(4:4:end),y4(4:4:end),'o','color','k','linewidth',2,'markersize',10,'markerfacecolor','w')
%% log scale along x:
%set(gca,'xsca','log')
%% set axes limits [xmin xmax ymin ymax]:
axis([0 2 -1.2 1.2])
%% x-axis ticks:
set(gca,'xtick',[0 1 2])
%% y-axis ticks:
set(gca,'ytick',[-1 0 1])
%% x and y labels:
xlabel('\gamma','position',[1.5 -1.35])
ylabel('\psi','position',[-.1 .5])
%%
box on
%% For a vectorial figure (for latex):
print -deps2c figure_name.eps
%% For an image figure at resolution 300 pixels/inch (for word):
%print -dpng -r300 figure_name.png
