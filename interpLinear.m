function [ y ] = interpLinear( xVec, yVec, x )
% Linear interpolation function. Faster than Matlab's interp1 with the
% 'linear' option for certain cases.
%
% Rick Fenrich 10/6/15

y = zeros(length(x),1);

for ii = 1:length(x)

    ind1 = find(xVec <= x(ii), 1, 'last');
    ind2 = find(xVec >= x(ii), 2, 'first');
    if(numel(ind1) == 0) % out of bounds
        x1 = x(ii);
        y1 = yVec(1);
    else % in bounds
        x1 = xVec(ind1(1));
        y1 = yVec(ind1(1));
    end
    
    if(numel(ind2) == 0) % out of bounds
        x2 = x(ii);
        y2 = yVec(end);
    else
        x2 = xVec(ind2(1));
        y2 = yVec(ind2(1));
    end    
    
    if(x1 == x2)
        y(ii) = y1;
    else
        y(ii) = (x(ii)-x1)*(y2-y1)/(x2-x1) + y1;
    end
    
end

end

