function s=quadeq(p,d,delta)
% Author: Marco Maffezzoli. Ver. 1.0.0, 11/2012.
%

a=d.'*d;
b=2*p.'*d;
c=p.'*p-delta^2;
z=(-b+sqrt(b^2-4*a*c))/(2*a);
s=p+z*d;

end