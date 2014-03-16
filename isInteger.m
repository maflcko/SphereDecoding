function [isInt] = isInteger(x)
%ISINTEGER( x ) if x is a scalar integer
%   Returns true if x is a scalar integer
isInt = ~ isempty(x) ...
    && isnumeric(x) ...
    && isreal(x) ...
    && length(size(x)) == 2 ...
    && all(size(x) == [1, 1]) ...
    && all(isfinite(x)) ...
    && all(x == fix(x));
end
