function u_closest = decode3(y, H)
%DECODE3( y, H ) Lattice Decode
%   u_closest = DECODE3( y, H ) returns the indexes of the closest lattice
%   vector for a given inverse lower triangular generator matrix H and any
%   row vector y.
%   (Algorithm 3 by A. Ghasemmehdi, E. Agrell)

% A. Ghasemmehdi, E. Agrell, "Faster Recursions in Sphere Decoding,"
% IEEE Trans. Inf. Theory, vol. 57, no. 6, pp. 3530-3536, June 2011.

validateInput(H, y, 0, 0);
n = size(H, 1);

%= Declaration =%
dist = zeros(1, n + 1);
u = zeros(1, n);
E = zeros(n, n);
step = zeros(1, n);
u_closest = zeros(1, n);

%= Initialization =%
rho_n = Inf;
k = n;
dist(n + 1) = 0;
E(n, :) = y * H;
u(n) = round(E(n, n));
gamma = (E(n, n) - u(n)) / H(n, n);
step(n) = sgn(gamma);
dist(n) = gamma ^ 2;
loop_down = true;

while (true)
 
    while (loop_down)
        if (~ (k == 1))
            k = k - 1;
            E(k, 1:k) = E(k + 1, 1:k) - gamma * H(k + 1, 1:k);
            u(k) = round(E(k, k));
            gamma = (E(k, k) - u(k)) / H(k, k);
            step(k) = sgn(gamma);
            dist(k) = dist(k + 1) + gamma ^ 2;
        else
            u_closest = u;
            rho_n = dist(1);
        end
        loop_down = (dist(k) < rho_n);
    end
 
    while (~ loop_down)
        if (k == n)
            return;
        else
            k = k + 1;
            u(k) = u(k) + step(k);
            step(k) = - step(k) - sgn(step(k));
            gamma = (E(k, k) - u(k)) / H(k, k);
            dist(k) = dist(k + 1) + gamma ^ 2;
        end
        loop_down = (dist(k) < rho_n);
    end
 
end
