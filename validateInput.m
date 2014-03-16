function validateInput(M, y, U_min, U_max)
%VALIDATEINPUT validates the input of the sphere decoding algorithms
%   Throws an error if the input can not be decoded using any of the
%   sphere decoding algorithms.

% Marco Falke [2014]

n = size(M, 1);

if (~ isInteger(U_min)) || (~ isInteger(U_max))
    error('decode:notInteger', ...
        'U_min and U_max have to be integer')
end

if (U_min > U_max)
    error('decode:inverseConstellation', ...
        'U_min can not be greater than U_max')
end

if any(diag(M) < 0)
    error('decode:negtiveDiag', ...
        'All diagonal elements must be positive')
end

if ~ isLowerTriangular(M)
    error('decode:notTriangular', ...
        'The input matrix has to be lower triangular')
end

if ~ (rank(M) == n)
    error('decode:linearDependent', ...
        'The input matrix has to be linearly independent')
end

if ~ (length(y) == n)
    error('decode:dimensionMismatch', ...
        'The input matrix and vector have to be of the same dimension')
end

end
