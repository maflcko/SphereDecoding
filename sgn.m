function signum = sgn(v)
signum = sign(sign(v)-.1); %return negative if v==0
return;