function [x,er,y]=trustsolve(fname,x,opt,varargin)
% TrustSolve - Numerical solution of nonlinear equations using
%              the Newton Trust Region with Dogleg or CG-Steihaug methods.
%
% [x,er,y]=trustsolve(fname,x,opt,z1,z2,...zn)
%
% Input:  fname - name of the m-file implementing f(x)
%         x - initial guess
%         opt - vector of options
%               opt(1): TolX in stopping rule (>eps)
%               opt(1): TolY in convergence test (>eps)
%               opt(3): Maximum number of iterations (>1)
%               opt(4): Max. radius of trust region (>10)
%               opt(5): Min. acceptable relative decrease in f(x) (>1e-4)
%               opt(6): Algorithm
%                       1: Dogleg
%                       2: CG-Steihaug
%               (if opt=[] or missing then the following default values
%                are used: TolX=sqrt(eps), TolY=sqrt(eps), MaxIer=500,
%                Delta=1e10, Eta=0.05, Alg=1)
%         z1,z2,...,zn - optional arguments passed to fname
%
% Stopping rules: (1) norm(y)<=TolY (convergence to a solution)
%                 (2) norm(step)<=TolX*(1+norm(x))
%                     (stagnation: the algorithm reports success iif
%                      norm(y)>TolY*(1+norm(y0)))
%                 (3) the maximum number of iterations has been reached
%
% Output: x - solution of f(x)=0
%         er - error code
%              0: OK, Convergence to a solution
%              1: Convergence to a point that is not a solution
%              2: Maximum number of iterations reached
%         y - value of f(x) at x
%
% Notes: please note that the algorithm will converge
%        (under some regularity coditions) to a critical point of
%        the funtion ||f(x)||^2, i.e. it will NOT necessarily converge to a
%        solution of f(x)!
%
% External calls: jacob1.m (one-sided finite diff. Jacobian)
%                 steihaug.m (CG solution of Trust region subproblem)
%                 dogleg.m (Dogleg solution of Trust region subproblem)
%
% Author: Marco Maffezzoli (marco.maffezzoli@unibocconi.it).
% Ver. 2.0.1, 11/2012.
% Based on the the algorithms described in Nocedal and Wright (1999),
% "Numerical Optimization", Springer, p. 299-300.
%
% This program is provided "as-is," without any express or implied
% warranty. Please, handle it with care!
%

if nargin<2
    error('I need at least two input arguments!')
elseif nargin==2||isempty(opt)
    tolx=sqrt(eps);
    toly=sqrt(eps);
    maxeval=500;
    dbar=1e10;
    eta=0.05;
    alg=1;
else
    tolx=max(opt(1),eps);
    toly=max(opt(2),eps);
    maxeval=max(round(opt(3)),1);
    dbar=max(opt(4),1);
    eta=max(opt(5),1e-4);
    alg=max(round(opt(6)),1);
end
shx=size(x);
% objf=str2func(fname);
objf=fname;
y=objf(x,varargin{:});
if numel(x)~=numel(y)
    error('The system is not square!')
end
ny=norm(y(:));
ny0=ny;
for k=1:maxeval
    J=jacob1(objf,x,y,varargin{:});
    H=J.'*J;
    g=J.'*y(:);
    if k==1
        ncp=norm(g)^3/(g.'*H*g);
        delta=max(1,min(ncp,dbar));
    end
    switch alg
        case 1
            [s,ns]=dogleg(g,H,delta);
        otherwise
            [s,ns]=steihaug(g,H,delta);
    end
    if ns<=tolx*(1+norm(x(:)));
        if ny<=toly*(1+ny0)
            er=0;
        else
            er=1;
        end
        return
    end
    x1=reshape(x(:)+s,shx);
    y1=objf(x1,varargin{:});
    ny1=norm(y1(:));
    rho=(ny1^2-ny^2)/(2*s.'*g+s.'*H*s);
    if rho<0.25
        delta=delta/4;
    elseif (rho>0.75)&&(ns>=delta)
        delta=min(2*delta,dbar);
    end
    if rho>eta
        x=x1;
        y=y1;
        ny=ny1;
    end
end
er=2;

end