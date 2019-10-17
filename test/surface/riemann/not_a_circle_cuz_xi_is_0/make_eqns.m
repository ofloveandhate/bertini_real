% a script for making equations for a circle as a reimann surface

syms xi xr yi yr
assume(xi,'real')
assume(xr,'real')
assume(yr,'real')
assume(yi,'real')

xi = 0;

x = xr+1i*xi;
y = yr+1i*yi;


f = x^2 + y^2 - 1;


f1 = real(f);
f2 = imag(f);


disp(sprintf('f1 = %s;',char(f1)))
disp(sprintf('f2 = %s;',char(f2)))