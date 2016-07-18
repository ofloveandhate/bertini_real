
function hbarplot=barplot(x,y,z,varargin)
% plot a 3d surface/tube centered at (x,y,z)
% usage:
%           hbarplot=barplot(x,y,z,varargin)
%           varargin:
%                       'v' : the cdata values, same size as z
%                       'n' : N-isogon as the basic shape
%                       'r' : radius of the N-isogon
%    e.g.,
%          zz=0:0.5:10;
%          hh=barplot(cos(zz/10*pi),5,zz,'v',randn(1,numel(zz))+sin(zz/10*pi),'n',4,'r',sin(zz/10*pi)); view(3)
%          http://scriptdemo.blogspot.com

%http://scriptdemo.blogspot.com/2012/08/maltab-bartube-plot-along-3d-path.html


myN=3; myR=0.5; isV=0;
if nargin>3
   while(size(varargin,2)>0)
      switch lower(varargin{1})
      case {'v','cdata','colordata','facecolordata'}
          v=varargin{2};
          varargin(1:2)=[];
          isV=1;
      case {'n','nshape'}
          myN=varargin{2};
          varargin(1:2)=[];
      case {'r','radius'}
          myR=varargin{2};
          varargin(1:2)=[];
      otherwise
          error(['Not defined properties!',varargin{1}])
      end
   end
elseif nargin==0
   % demo case
   zz=0:0.5:10;
   hbarplot=barplot(cos(zz/10*pi),5,zz,'v',randn(1,numel(zz))+sin(zz/10*pi),'n',4,'r',sin(zz/10*pi));
   view(-150,30);
   return
elseif nargin<3
   help barplot
   return
end

[xx,yy,zz]=getMeshN(x,y,z,myR,myN);
if isV==0
   hbarplot=surface(xx',yy',zz');
else
   if numel(v)==numel(z)
      vv=repmat(reshape(v,[],1),1,size(xx,2));
   end
   hbarplot=surface(xx',yy',zz');
   set(hbarplot,'cdata',vv');
end
set(hbarplot,'tag','barplot','linestyle','none','facecolor','interp');

set(gcf,'color','w');






function [xx,yy,zz]=getMeshN(x0,y0,z,myR,N)
if nargin==3
   myR=0.5;
   N=3;
elseif nargin==4
   N=3;
elseif nargin~=5
   error('Usage: [xx,yy,zz]=getMeshN(x,y,z,myR,N)');
end
z=reshape(z,[],1);
numZ=numel(z);

myAng=linspace(0,2*pi,N+1);
xslab=cos(myAng); yslab=sin(myAng);
%N=numel(xslab)-1; % could be other shape as well.

if numel(x0)==1
   x0=repmat(x0,numZ,N+1);
elseif numel(x0)==numZ
   x0=repmat(reshape(x0,[],1),1,N+1);
end
if numel(y0)==1
   y0=repmat(y0,numZ,N+1);
elseif numel(y0)==numZ
   y0=repmat(reshape(y0,[],1),1,N+1);
end

if numel(myR)==1
   xx=repmat(xslab*myR,numZ,1)+x0;
   yy=repmat(yslab*myR,numZ,1)+y0;
else
   myR=repmat(reshape(myR,[],1),1,N+1);
   xx=repmat(xslab,numZ,1).*myR+x0;
   yy=repmat(yslab,numZ,1).*myR+y0;
end
zz=repmat(z,1,N+1);
