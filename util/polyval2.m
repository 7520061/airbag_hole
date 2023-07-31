function res = polyval2(coef,x)
    order = length(coef)-1;
    exponent = repmat([order:-1:0]', 1, length(x));
    res = coef * x.^exponent;
end