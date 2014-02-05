function u_closest = decode8( y, H , U_min, U_max)
%DECODE8( y, H ) Lattice Decode
%   u_closest = DECODE8( y, H ) returns the indexes of the closest lattice
%   vector for a given inverse lower triangular generator matrix H and any
%   row vector y.
%   (Algorithm 8 by A. Ghasemmehdi, E. Agrell)

% A. Ghasemmehdi, E. Agrell, "Faster Recursions in Sphere Decoding,"
% IEEE Trans. Inf. Theory, vol. 57, no. 6, pp. 3530-3536, June 2011.

validateInput(H,y,U_min,U_max);

%= Initialize =%
n=size(H,1);
rho_n=Inf;
k=n+1;
d=zeros(1,n)+n;
dist(n+1)=0;
E(n,:)=y*H;
loop_down=true;

gamma(n)=0; %init needed due to matlab 'restrictions'

while(true) %LOOP_LABEL
    
    while(loop_down)
        if(~(k==1))
            k=k-1;
            for a=d(k):-1:(k+1)
                E(a-1,k)=E(a,k)-gamma(a)*H(a,k);
            end
            u(k)=roundc(E(k,k),U_min,U_max);
            gamma(k)=(E(k,k)-u(k))/H(k,k);
            step(k)=sgn(gamma(k));
            dist(k)=dist(k+1)+gamma(k)^2;
        else
            u_closest=u;
            rho_n=dist(1);
        end
        loop_down=(dist(k)<rho_n);
    end
    m=k;
    while(~loop_down)
        if (k==n)
            return;
        else
            k=k+1;
            gamma(k)=Inf;
            u(k)=u(k)+step(k);
            step(k)=-step(k)-sgn(step(k));
            if ( (U_min<=u(k)) && (u(k)<=U_max) )
                gamma(k)=(E(k,k)-u(k))/H(k,k);
            else
                u(k)=u(k)+step(k);
                step(k)=-step(k)-sgn(step(k));
                if ( (U_min<=u(k)) && (u(k)<=U_max) )
                    gamma(k)=(E(k,k)-u(k))/H(k,k);
                end
            end
            dist(k)=dist(k+1)+gamma(k)^2;
        end
        loop_down=(dist(k)<rho_n);
    end
    
    d(m:k-1)=k;
    for a=(m-1):-1:1
        if (d(a)<k)
            d(a)=k;
        else
            break;
        end
    end
    
end % GOTO LOOP_LABEL
