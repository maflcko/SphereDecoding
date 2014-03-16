function [s] = roundc(t, U_min, U_max)
%ROUNDC used in finite constellations
%   S = ROUNDC( t, U_min, U_max ) returns the index of the nearest elegible
%   layer

s = min(U_max, max(U_min, round(t)));

end
