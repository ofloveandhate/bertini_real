function saveobjtube(fname,x,y,z,varargin)

% Save a tube plot of a number of tubes as a Wavefront .obj file 
%
% saveobjtube(filename, x, y, z) - Saves a tube plot of the tubes
% in the NxT vectors.
%
% saveobjtube(filename,x,y,z,radius,param,subdivs,vec) - As tubeplot (see function)
% x,y,z denotes the path, radius its radius, param parametrisation 
% (see below), subdivs the number of subdivisions and vec an optional 
% alternate frame.
%
% Different parametrisations
%
% Type 1: each tube parametrised along length and circumference to
% unit square
%
% Type 2: each tube parametrised along length and tube # to unit square
%
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005

N=size(x,1);
T=size(x,2);

param=1;
subdivs = 6;
r=0.1;

if (nargin < 5)
  r=x*0+1;
else
  r=varargin{1};
end
if (nargin > 5)
  param=varargin{2};
end
if (nargin > 6)
  subdivs=varargin{3}+1;
end
if (nargin > 7)
  vec=varargin{4};
end

fid=fopen(fname,'w');

nn=1;

for t=1:T  
  fprintf(fid,'g tube%d\n',t);

  if (T>1)
    xx=x(:,t);
    yy=y(:,t);
    zz=z(:,t);
    rr=r(:,t);
  else
    xx=x;
    yy=y;
    zz=z;
    rr=r;
  end
  
  if (nargin > 7)
    vec=varargin{4};
    [X,Y,Z]=tubeplot(xx,yy,zz,rr,rr,subdivs,vec);
  else
    [X,Y,Z]=tubeplot(xx,yy,zz,rr,rr,subdivs);
  end
  
  l=size(X,1); h=size(X,2);  
  n=zeros(l,h);
  for i=1:l
    for j=1:h
      n(i,j)=nn; 
      fprintf(fid, 'v %f %f %f\n',X(i,j),Y(i,j),Z(i,j)); 
      if (param==1)
	fprintf(fid, 'vt %f %f\n',(i-1)/(l-1),(j-1)/(h-1)); 
      else
	fprintf(fid, 'vt %f %f\n',(i-1)/(l-1),(t-1)/(T)+((j-1)/h)/T); 
      end
      nn=nn+1;
    end
  end
  
  for i=1:(l-1)
    for j=1:(h-1)
      fprintf(fid,'f %d/%d %d/%d %d/%d %d/%d\n',n(i,j),n(i,j),n(i+1,j),n(i+1,j),n(i+1,j+1),n(i+1,j+1),n(i,j+1),n(i,j+1));
    end
  end
  fprintf(fid,'g\n');
end

fclose(fid);
