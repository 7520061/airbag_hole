function [Yq] = calc1DInterp(X, Y, Xq, pp)

Yq = lininterp1f(X, Y, Xq, 1)';   

% Slower alternatives
% Yq = pp(Xq);    
% Yq = interp1(X, Y, Xq);    
% Yq = interp1qr(X, Y, Xq);   

end