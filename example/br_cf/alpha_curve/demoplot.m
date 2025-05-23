function C = demoplot(C, f)

if nargin==0
	C = BrCfCurve();
	C.ComputeChebfuns();
end
if nargin <= 1
	f = @(x1,x2) besselj(1,(abs(x1+x2))+(erfc(x1-x2)-1));
end
b = 18;
x = linspace(-b,b,1000);y = x;
[X,Y] = meshgrid(x,y);

figure(1)
pcolor(X,Y,f(X,Y)); shading interp; hold on

curve_handles = plot(C);
set(curve_handles,'LineWidth',8,'Color','w');

curve_handles = plot(C);
set(curve_handles,'LineWidth',2,'Color','k');

C.set_function(f);
[max_vals,max_locations] = max(C,'local');

max_handles = plot(max_locations(:,1),max_locations(:,2));
set(max_handles,'MarkerSize',14,'LineStyle','none','Marker','o','MarkerFaceColor',[0.94 0.94 0.94],'LineWidth',2,'Color','k');

axis square;


set(gca,'FontSize',14);
xlabel('x');
ylabel('y');




figure(2)

pcolor(X,Y,f(X,Y)); shading interp; hold on

curve_handles = plot(C);
set(curve_handles,'LineWidth',8,'Color','w');

curve_handles = plot(C);
set(curve_handles,'LineWidth',2,'Color','k');


[max_vals,max_locations] = max(C);

max_handles = plot(max_locations(:,1),max_locations(:,2));
set(max_handles,'MarkerSize',14,'LineStyle','none','Marker','o','MarkerFaceColor',[0.94 0.94 0.94],'LineWidth',2,'Color','k');

axis([-b, b, -b, b])
axis square;


set(gca,'FontSize',14);
xlabel('x');
ylabel('y');

end