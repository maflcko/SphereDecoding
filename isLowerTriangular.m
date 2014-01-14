function out = isLowerTriangular( M )
%ISLOWERTRIANGULAR True for lower triangular matrix
%   ISLOWERTRIANGULAR( M ) returns false if the matrix M has non-zero 
%   entries in the triangle above the diagonal, and true otherwise

% Marco Falke [2014]

out=isUpperTriangular(M');
return;

end

