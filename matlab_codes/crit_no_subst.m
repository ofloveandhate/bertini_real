% [] = crit_no_subst(filename,pi_vals,OutputName)
%
% computes a critical curve from an input system.
%
% filename - the name of the file to parse and deflate
% pi_vals - 
% OutputName - the name of the file to create and write to
%
% dani brake
% 2017

function crit_no_subst(filename,pi_vals,OutputName)

b_input = bertini_input(filename); % from the Bertini_tropical package

system_variables = b_input.variable_group;
num_vars = length(system_variables);

num_funcs = size(b_input.functions,1);

vars = sym(zeros(1,num_vars));
for ii = 1:num_vars
	eval(sprintf('syms %s',system_variables{ii}));  %make the variable be a symbol in matlab memory
	vars(ii) = sym(system_variables{ii});
	%vars is used in the jacobian call, so that only those derivatives are computed
end

if ~isempty(b_input.constant)
    system_constants = b_input.constant{:,1};
    num_constants = length(system_constants);

    for ii = 1:num_constants
        eval(sprintf('syms %s',system_constants(ii)));
    end
end

ensure_fns_are_subfns(b_input);


if ~isempty(b_input.subfunction)
	subfunc_names = b_input.subfunction(:,1);
else
	subfunc_names = {};
end
num_subfuncs = length(subfunc_names);


if num_subfuncs > 0
	[var_deps,subfunc_deps] = dependency_graph(b_input); % determine which subfunctions depend on what, computed by regular expressions
end


absolute_var_deps = zeros(num_subfuncs, num_vars);
for ii= 1:num_subfuncs
	for jj = 1:num_vars
		absolute_var_deps(ii,jj) = depends_on_var(ii,jj, var_deps, subfunc_deps);
	end
end



orig_subfunc_args = cell(num_subfuncs,2);

fprintf('\tsetting up subfunctions in memory\n')

for ii = 1:num_subfuncs
	curr_subfunc_args = '';
	curr_subfunc_args_for_regexp = '';
	for jj = 1:num_vars
		if absolute_var_deps(ii,jj)
			curr_subfunc_args = sprintf('%s%s, ',curr_subfunc_args,system_variables{jj});
			curr_subfunc_args_for_regexp = sprintf('%s%s,\\s*',curr_subfunc_args_for_regexp,system_variables{jj});
		end
	end
	curr_subfunc_args = curr_subfunc_args(1:end-2);
	curr_subfunc_args_for_regexp = [curr_subfunc_args_for_regexp(1:end-4) '\s*'];


	orig_subfunc_args{ii,1} = curr_subfunc_args; %store it.
	orig_subfunc_args{ii,2} = curr_subfunc_args_for_regexp; %store it.

	currstr = sprintf('syms %s(%s)',b_input.subfunction{ii,1},curr_subfunc_args);
	eval(currstr);
end



%this block creates the functions, preserving subfunctions as subfunctions,
%without substitution.

fprintf('\tmaking functions in memory\n');
function_names = b_input.functions(:,1);
f_user = sym(zeros(length(function_names),1));
for ii = 1:length(function_names)
	curr_func = b_input.functions{ii,2};
	for jj = 1:num_subfuncs
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{jj});
		newname = sprintf('%s(%s)',subfunc_names{jj},orig_subfunc_args{jj,1});
		newpattern = sprintf('$1%s$2',newname);
		curr_func = regexprep(curr_func,oldpattern,newpattern);
	end
	f_user(ii) = eval(curr_func);
end

Jac = jacobian(f_user,vars);


% initialize to empty
new_subfunc_names = {};
is_zero_derivative = zeros(length(subfunc_names),num_vars);

fprintf('\tcomputing partial derivatives of subfunctions\n');
for ii = 1:length(subfunc_names)
	for jj = 1:num_vars

		if absolute_var_deps(ii,jj)~=0

			curr_subfunc = b_input.symbol_value(subfunc_names{ii});

			for kk = 1:num_subfuncs
				oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{kk});
				newname = sprintf('%s(%s)',subfunc_names{kk},orig_subfunc_args{kk,1});
				newpattern = sprintf('$1%s$2',newname);
				curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
			end

			current_derivative = diff(eval(curr_subfunc),system_variables{jj});
			if current_derivative==0
				is_zero_derivative(ii,jj) = 1;
			end
			symname = sprintf('DIFF_%s_%s',subfunc_names{ii},system_variables{jj});
			if ~b_input.is_symbol_declared(symname)
				b_input.declare_and_define(symname,current_derivative,'subfunction');
				new_subfunc_names{end+1} = symname;
			end
		else
			is_zero_derivative(ii,jj) = 1;
		end


	end
end

%add two rows of pi symbols onto jacobian, take determinant (through 161)
num_pis = num_vars - num_funcs;

%bring in pi values and fill Jacobian out with them
for ii = 2:num_pis+1
    for jj = 1:num_vars
        
        new_pi = sprintf('pi_%i_%i',ii,jj);
        new_val = pi_vals{ii-1,jj};
        
        b_input.declare_and_define(new_pi,new_val,'constant');
        
        Jac(ii,jj) = new_pi;
    end
end

crit_curve_fn = sym([]);

temp_fn = det(Jac);
if temp_fn ~= 0
    crit_curve_fn = temp_fn;
end



%transform from loop to single statement
fprintf('\tdoing partial derivative substitutions for detjac\n')

curr_crit_curve = char(crit_curve_fn);
    base = sprintf('diff\\(%s\\(%s\\)',subfunc_names{1},orig_subfunc_args{1,2});
    for kk = 1:num_vars

        oldname = sprintf('%s,\\s*%s\\)',base,system_variables{kk});
        oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
        newname = sprintf('DIFF_%s_%s',subfunc_names{1},system_variables{kk});
        newpattern = sprintf('$1%s$2',newname);
        curr_crit_curve = regexprep(curr_crit_curve,oldpattern,newpattern);
    end
crit_curve_fn = curr_crit_curve;




fprintf('\tdoing subfuncion substitutions for defl_fns\n')

%now we substitute away the subfunction f(...) statements, by the names of
%the original subfunctions.
%transform from loop to single statement

curr_crit_curve = char(crit_curve_fn);
    oldname = sprintf('%s\\(%s\\)',subfunc_names{1},orig_subfunc_args{1,2});
    oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
    newname = sprintf('%s',subfunc_names{1});
    newpattern = sprintf('$1%s$2',newname);
    curr_crit_curve = regexprep(curr_crit_curve,oldpattern,newpattern);
    
crit_curve_fn = sprintf('(%s)',curr_crit_curve);




% finally, we need to substitute away any remaining partial derivative
% statements in the subfunctions.
orig_subfunc_names = subfunc_names;

fprintf('\tdoing substitutions for new subfunctions\n');

for ii = length(new_subfunc_names):-1:1

	ind = b_input.symbol_index(new_subfunc_names{ii},'subfunction');
	curr_subfunc = char(b_input.subfunction{ind,2});
	b_input.subfunction{ind,2} = sub_away_diffs(curr_subfunc, orig_subfunc_names,orig_subfunc_args, system_variables);
end

crit_curve_fn = sub_away_diffs(crit_curve_fn, orig_subfunc_names, orig_subfunc_args, system_variables);

b_input.declare_and_define(sprintf('crit_curve'),crit_curve_fn,'function');


if nargin == 4
    output_filename = sprintf('%s_deflated_%i',filename,defl_iteration);
else
    output_filename = OutputName;
end

write_bertini_input_file(b_input.variable_group, b_input.functions,'filename',output_filename,'options',b_input.config,'constants',b_input.constant,'subfunctions',b_input.subfunction);


end



function [curr_subfunc] = sub_away_diffs(curr_subfunc, orig_subfunc_names, orig_subfunc_args, system_variables)

num_vars = length(system_variables);

    for jj = length(orig_subfunc_names):-1:1
		for kk = 1:num_vars
			oldname = sprintf('diff\\(%s\\(%s\\),\\s*%s)',orig_subfunc_names{jj},orig_subfunc_args{jj,2},system_variables{kk});
			oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);

			newname = sprintf('DIFF_%s_%s',orig_subfunc_names{jj},system_variables{kk});
			newpattern = sprintf('$1%s$2',newname);
			curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
		end

		oldname = sprintf('%s\\(%s\\)',orig_subfunc_names{jj},orig_subfunc_args{jj,2});
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);

		newname = orig_subfunc_names{jj};
		newpattern = sprintf('$1%s$2',newname);
		curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
    end
    
end

function ensure_fns_are_subfns(b_input)

for ii = 1:size(b_input.functions,1)
    
    is_subf=false;
    fn_val = b_input.functions{ii,2};
    fn_name = b_input.functions{ii,1};
    
    for jj = 1:size(b_input.subfunction,1)
        
        subfn_name = b_input.subfunction{jj,1};
        
        if strcmp(strtrim(fn_val),strtrim(subfn_name))
            is_subf=true;
            break;
            
        end
    end
    
    if is_subf == true
        
    elseif is_subf == false
        
        new_subf_sym = sprintf('subf_%s',fn_name);
        f_sym = sprintf('%s',fn_name);
        sym_val = b_input.functions{ii,2};
        
        b_input.declare_and_define(new_subf_sym,sym_val,'subfunction');
        b_input.define_symbol(f_sym,new_subf_sym);
    end
end

end