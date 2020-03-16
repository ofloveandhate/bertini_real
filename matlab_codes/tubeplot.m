function [varargout]=tubeplot(x,y,z,varargin)  

% TUBEPLOT - plots a tube r along the space curve x,y,z.
%
% tubeplot(x,y,z) plots the basic tube with radius 1
% tubeplot(x,y,z,r) plots the basic tube with variable radius r (either a vector or a value)
% tubeplot(x,y,z,r,v) plots the basic tube with coloring dependent on the values in the vector v
% tubeplot(x,y,z,r,v,s) plots the tube with s tangential subdivisions (default is 6)
%
% [X,Y,Z]=tubeplot(x,y,z) returns [Nx3] matrices suitable for mesh or surf
%
% Note that the tube may pinch at points where the normal and binormal 
% misbehaves. It is suitable for general space curves, not ones that 
% contain straight sections. Normally the tube is calculated using the
% Frenet frame, making the tube minimally twisted except at inflexion points.
%
% To deal with this problem there is an alternative frame:
% tubeplot(x,y,z,r,v,s,vec) calculates the tube by setting the normal to
% the cross product of the tangent and the vector vec. If it is chosen so 
% that it is always far from the tangent vector the frame will not twist unduly
%
% Example:
%
%  t=0:(2*pi/100):(2*pi);
%  x=cos(t*2).*(2+sin(t*3)*.3);
%  y=sin(t*2).*(2+sin(t*3)*.3);
%  z=cos(t*3)*.3;
%  tubeplot(x,y,z,0.14*sin(t*5)+.29,t,10)
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005
% 
% acquired from
% http://www.aleph.se/Nada/Ray/Tubeplot/tubeplot.html
% 
% Danielle sent many many emails to many people, looking 
% for explicit permission from the author to mod & distribute this code.
% There is no license attached, but the website used for distribution
% says "The package is free for anybody to use."  
% I think they meant free as in free speech, not free beer.
% I did my due diligence to acquire permission, ultimately interpreting the
% online statement to permit me to redistribute it.  

	opt = parse_args(varargin{:});

  n = opt.n;

  N=size(x,1);
  if (N==1)
    x=x';
    y=y';
    z=z';
    N=size(x,1);
  end
	
  r=x*0+opt.radius;
  

   [t_frame,n_frame,b_frame]=frenet(x,y,z);


  

  
  vertices = zeros(0,3);



	theta=linspace(0,2*pi,n+1);
	theta=theta(1:end-1); % cuz its periodic, yo
  for i=1:N
    x_=x(i) + r(i)*(n_frame(i,1)*cos(theta) + b_frame(i,1)*sin(theta))';
    y_=y(i) + r(i)*(n_frame(i,2)*cos(theta) + b_frame(i,2)*sin(theta))';
    z_=z(i) + r(i)*(n_frame(i,3)*cos(theta) + b_frame(i,3)*sin(theta))';
	
	vertices = [vertices;[x_ y_ z_]];

  end

  faces = [];
  for i=1:N-1
	  
	ell = (i-1)*(n)+1:( (i)*(n));
	
	a = ell;
	b = circshift(a,[0 1]);
	c = a+n;
	
	t1 = [b;...
		a;...
		c];
	
	t2 = [b;...
		c;...
		circshift(c,[0 1])];
	faces = [faces;t1';t2'];
  end

  fv.vertices = vertices;
  fv.faces = faces;
  
if opt.closed_left
	fv.vertices = [fv.vertices;[x(1) y(1) z(1)]];
	t1 = [[2:n 1]; 1:n; ones(1,n)*length(fv.vertices(:,1))];
	fv.faces = [fv.faces; t1'];
end


if opt.closed_right
	fv.vertices = [fv.vertices;[x(end) y(end) z(end)]];
	t1 = [(1:n)+(N-1)*(n); [(2:n) 1]+(N-1)*(n); ones(1,n)*length(fv.vertices(:,1))];
	fv.faces = [fv.faces; t1'];
end


if opt.render
	h = patch(fv);
	set(h,'FaceAlpha',0.2);
	set(h,'FaceColor',[1 0 0]);
	set(h,'EdgeColor',[0.5 0 0]);
end


fv.left_cap = 1:n;
fv.right_cap = (1:n)+(N-1)*(opt.n);

varargout{1} = fv; 
    
  end

  
function opt = parse_args(varargin)
		
opt.radius = 0.15;
opt.n = 100;

opt.render_curves = true;
opt.render = true;
opt.write_to_stl = false;

opt.closed_left = false;
opt.closed_right = false;

while ~isempty(varargin)
	curr_opt = varargin{1};
	varargin = varargin(2:end);
	switch curr_opt
		case {'radius','r'}
			opt.radius = varargin{1};
			varargin = varargin(2:end);
		case 'n'
			opt.n = varargin{1};
			varargin = varargin(2:end);
		case 'closed'
			curr_val = varargin{1};
			varargin = varargin(2:end);
			switch curr_val
				case 'left'
					opt.closed_left = true;
				case 'right'
					opt.closed_right = true;
				case 'both'
					opt.closed_left = true;
					opt.closed_right = true;
				case 'none'
					opt.closed_left = false;
					opt.closed_right = false;
				otherwise
					error('bad option %s to option ''closed'' for tubeplot', curr_val);
			end
		case 'render'
			opt.render = varargin{1};
			varargin = varargin(2:end);
		otherwise
			error('bad option ''%s'' to tubeplot', curr_opt);
	end
end

end