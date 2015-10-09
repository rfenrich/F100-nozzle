function [s,ns]=steihaug(g,H,delta)
% Author: Marco Maffezzoli. Ver. 2.0.3, 05/2009.
%
% This function is strictly based on the algorithm described in Steihaug,
% "The Conjugate Gradient Method and Trust Regions in Large Scale
% Optimization," SIAM Journal of Numerical Analysis, 20(3), 1983. 
% See the article for details.
%

r=-g;
d=r;
ng=norm(g);
z=ng^2;
s=zeros(length(r),1);
tol=sqrt(eps)*ng;
ns=0;
while sqrt(z)>tol
    gamma=d.'*H*d;
    if gamma<=0
        s=quadeq(s,d,delta);
        ns=delta;
        return
    end
    alpha=z/gamma;
    s1=s+alpha*d;
    ns=norm(s1);
    if ns>=delta
        s=quadeq(s,d,delta);
        ns=delta;
        return
    end
    s=s1;
    r=r-alpha*H*d;
    z1=r.'*r;
    d=r+(z1/z)*d;
    z=z1;
end

end