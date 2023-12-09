% generates an evaluable or callable object, which produces the jacobian of the
% system.  huzzah.
%
% depends, like so many other bits of code i write, on Brakelab.
%
% works only for one variable group systems
%
% this code is admittedly *dumb* and *fragile*
%  
% i really do apologize for this code.
%
% silviana amethyst, 2019


function jac_evaluable = generate_evaluable_jacobian(BRinfo)
b = bertini_input(); % make an empty input file.  golly i am glad i have this written already.
b.parse_from_string(BRinfo.input);


num_vars = length(b.variable_group);
num_funcs = size(b.functions,1);


syms_str = 'syms';
for ii = 1:num_vars
	syms_str = sprintf('%s %s',syms_str,b.variable_group{ii});
end

% now evaluate the constants so they're in memory.
for ii = 1:size(b.constant,1)
	eval(sprintf('%s = %s;',b.constant{ii,1},b.constant{ii,2}));
end

eval(syms_str); % groovy, now we have some symbols in memory.

sys = sym([]);
for ii = 1:num_funcs
	func_str = sprintf('f = %s;',b.functions{ii,2});
	eval(func_str);
	sys(end+1) = f;
end

J = jacobian(sys);
J_as_char = cell(size(J));
[m,n] = size(J);
for ii = 1:m
	for jj = 1:n
		J_as_char{ii,jj} = replace_vars(J(ii,jj),b);
		
	end
end



J_handle_str = sprintf('jac_evaluable = @(x) [');
for ii = 1:m
	J_handle_str = sprintf('%s',J_handle_str);
	for jj = 1:n
		if jj ==n
			terminator = '';
		else
			terminator = ',';
		end
		J_handle_str = sprintf('%s %s%s',J_handle_str,J_as_char{ii,jj},terminator);
	end
	J_handle_str = sprintf('%s;',J_handle_str);
end
J_handle_str = sprintf('%s];',J_handle_str);

	
eval(J_handle_str); % makes the returned value
end

function J_str = replace_vars(J_entry,b)
J_str = char(J_entry);

num_vars = length(b.variable_group);
for ii = 1:num_vars
	pattern = sprintf('%s',b.variable_group{ii});
	replace = sprintf('x(%i)',ii);
	J_str = regexprep(J_str,pattern,replace);
end
end

