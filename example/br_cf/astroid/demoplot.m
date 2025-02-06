function C = demoplot(C, f)

if nargin==0
	C = BrCfCurve();
	C.ComputeChebfuns();
end
if nargin <= 1
	f = @(x1,x2) sin(7*x1)+sin(5*x2);
end
b = 1.5;
x = linspace(-b,b,1000);y = x;
[X,Y] = meshgrid(x,y);

pcolor(X,Y,f(X,Y)); shading interp; hold on

curve_handles = plot(C);
set(curve_handles,'LineWidth',8,'Color','w');

curve_handles = plot(C);
set(curve_handles,'LineWidth',2,'Color','k');


[max_vals,max_locations] = min(C,f,'local');

max_handles = plot(max_locations(:,1),max_locations(:,2));
set(max_handles,'MarkerSize',14,'LineStyle','none','Marker','o','MarkerFaceColor',[0.94 0.94 0.94],'LineWidth',2,'Color','k');

axis square;

xlabel('x');
ylabel('y');
set(gca,'FontSize',14);
end