function [y, G] = preprocess(x, D)
%PREPROCESS alorithm for lattice decode
%   [y, G] = PREPROCESS(x, D)
%   Returns the lower triangular matrix G and the vector y used for the
%   lattice decoding algorithms by Agrell et. al. Input is any generator
%   matrix D and the vector x which is to be decoded.
%   Hint: 
%     It is noted that the generator matrix must contain the basis
%     vectors in it's rows. Consequently x is a row-vector as well.
%     (Transpose your input otherwise)

[Q2, R2] = qr(D');
S=diag(sign(diag(R2)));
G=R2'*S;
Q=Q2*S;
y=x*Q;


return;
