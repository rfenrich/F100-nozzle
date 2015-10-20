function [ y ] = interpLinear( xVec, yVec, x )
% Linear interpolation function. Faster than Matlab's interp1 with the
% 'linear' option for certain cases.
%
% Rick Fenrich 10/6/15

y = zeros(length(x),1);

for ii = 1:length(x)

    ind1 = find(xVec <= x(ii), 1, 'last');
    ind2 = find(xVec >= x(ii), 2, 'first');
    x1 = xVec(ind1(1));
    x2 = xVec(ind2(1));
    y1 = yVec(ind1(1));
    
    if(x1 == x2)
        y(ii) = y1;
    else
        y2 = yVec(ind2(1));
        y(ii) = (x(ii)-x1)*(y2-y1)/(x2-x1) + y1;
    end
    
end

end

