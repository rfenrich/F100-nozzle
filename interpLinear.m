function [ y ] = interpLinear( xVec, yVec, x )
% Linear interpolation function
%
% Rick Fenrich 10/6/15

y = zeros(length(x));

for ii = 1:length(x)

    ind1 = find(xVec <= x(ii), 1, 'last');
    ind2 = find(xVec >= x(ii), 2, 'first');
    x1 = xVec(ind1(1));
    x2 = xVec(ind2(1));
    y1 = yVec(ind1(1));
    y2 = yVec(ind2(1));

    y(ii) = (x(ii)-x1)*(y2-y1)/(x2-x1) + y1;
    
end

end

