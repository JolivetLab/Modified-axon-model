function Y = limitCases(X, Y, r, T)

I = isinf(Y)|isnan(Y);
if ~any(I), return, end
I = find(I);
for i = 1 : length(I)
    x = X(I(i))-1e-6:1e-7:X(I(i))+1e-6;
    y = r(x, T);
    Y(I(i)) = interp1(x(x~=X(I(i))), y(x~=X(I(i))), X(I(i)));
end