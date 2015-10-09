function J=jacob1(objf,x,y,varargin)
% Jacob1 - One-sided finite difference Jacobian
%
% J=Jacob1(fname,x,y,z1,z2,...zn)
%
% Input:  objf - handle to the function implementing f(x)
%         x - point where to calculate the Jacobian
%         y - value of f(x) at x
%         z1,z2,...,zn - optional arguments passed to fname
%
%         NB: x and y can be vectors or matrices of ANY size
%
% Output: J - Jacobian of f(x) at x
%
%         NB: the Jacobian is a nxm matrix, where n is numel(x) and m is
%         numel(y). If g_k equals f(x+h*e_k), then g_k(:) is the kth columnt of J.
%
% Author: Marco Maffezzoli. Ver. 2.0.9, 11/2012.
%
% This program is provided "as-is," without any express or implied
% warranty. Please, handle it with care!
%

shx=size(x);
m=numel(y);
n=numel(x);
h=spdiags(sqrt(eps)*max(abs(x(:)),1),0,n,n);
J=zeros(m,n);
for k=1:n
    y1=objf(reshape(x(:)+h(:,k),shx),varargin{:});
    J(:,k)=(y1(:)-y(:))/h(k,k);
end

end