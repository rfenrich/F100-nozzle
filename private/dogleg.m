function [s,ns,err]=dogleg(g,H,delta)
% Author: Marco Maffezzoli. Version: 2.0.4, 11/2012.
%
% Compute the full Newton step if H is well conditioned, otherwise perform
% the perturbation suggested by Dennis and Schnabel, "Numerical Methods for
% Unconstrained Optimization and Nonlinear Equations," 1996, p. 151.
%
% Note the so-far undocumented matlab feature that temporarily converts a
% warning into an error, in order to allow the use of try-catch
%

w=warning('error','MATLAB:singularMatrix'); % undocumented feature!
try
    s=-H\g;
catch err
    n=length(g);
    s=-(H+realsqrt(n*eps)*norm(H,1)*eye(n))\g;
end
warning(w);
ns=norm(s);
if ns>delta
    ns=delta;
    ng=norm(g);
    ps=-delta*(g/ng);
    tau=min(1,ng^3/(delta*g.'*H*g));
    pc=tau*ps;
    if tau==1
        s=pc;
    else
        d=s-pc;
        s=quadeq(pc,d,delta);
    end
end

end