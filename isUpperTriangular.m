function out = isUpperTriangular( M )
%ISUPPERTRIANGULAR True for upper triangular matrix
%   ISUPPERTRIANGULAR( M ) returns false if the matrix M has non-zero
%   entries in the triangle below the diagonal, and true otherwise

% Marco Falke [2014]

out=false;
if~(size(M,1)==size(M,2))
    return;
end
maxAbs = max(max(abs(M)));
for row=2:size(M,1)
    if any(abs(M(row,1:row-1))>maxAbs*10^-12)
        return;
    end
end
out =true;

end

