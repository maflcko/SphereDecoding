function signum = sgn(v)
%SGN( v ) returns the signum of v
%   S = SGN( v ) returns the signum of v and negative if (v == 0)

signum = sign(sign(v) - .1);

end
