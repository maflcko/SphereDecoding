function u_closest = decode4( y, H , U_min, U_max)
%DECODE4( y, H ) Lattice Decode
%   u_closest = DECODE4( y, H ) returns the indexes of the closest lattice
%   vector for a given inverse lower triangular generator matrix H and any
%   row vector y.
%   (Algorithm 4 by A. Ghasemmehdi, E. Agrell)

% A. Ghasemmehdi, E. Agrell, "Faster Recursions in Sphere Decoding,"
% IEEE Trans. Inf. Theory, vol. 57, no. 6, pp. 3530-3536, June 2011.

n=size(H,1);

if any(diag(H)<0)
    error('decode:negtiveDiag',...
        'All diagonal elements must be positive')
end
if ~isLowerTriangular(H) 
    error('decode:notTriangular',...
        'The input matrix has to be lower triangular')
end
if ~(rank(H)==n)
    error('decode:linearDependent',...
        'The input matrix has to be linearly independent')
end
if ~(length(y)==n)
    error('decode:dimensionMismatch',...
        'The input matrix and vector have to be of the same dimension')
end

%= Initialize =%
rho_n=Inf;
k=n;
dist(n+1)=0;
E(n,:)=y*H;
u(n)=roundc(E(n,n),U_min,U_max);
gamma=(E(n,n)-u(n))/H(n,n);
step(n)=sgn(gamma);
dist(n)=gamma^2;

while(true) %LOOP_LABEL
    
    loop_down=true;
    while(loop_down)
        if(~(k==1))
            k=k-1;
            for a=1:k
                E(k,a)=E(k+1,a)-gamma*H(k+1,a);
            end
            u(k)=roundc(E(k,k),U_min,U_max);
            gamma=(E(k,k)-u(k))/H(k,k);
            step(k)=sgn(gamma);
            dist(k)=dist(k+1)+gamma^2;
        else
            u_closest=u;
            rho_n=dist(1);
        end
        loop_down=(dist(k)<rho_n);
    end
    
    loop_up=true;
    while(loop_up)
        if (k==n)
            return;
        else
            k=k+1;
            gamma=Inf;
            u(k)=u(k)+step(k);
            step(k)=-step(k)-sgn(step(k));
            if ( (U_min<=u(k)) && (u(k)<=U_max) )
                gamma=(E(k,k)-u(k))/H(k,k);
            else
                u(k)=u(k)+step(k);
                step(k)=-step(k)-sgn(step(k));
                if ( (U_min<=u(k)) && (u(k)<=U_max) )
                    gamma=(E(k,k)-u(k))/H(k,k);
                end
            end
            dist(k)=dist(k+1)+gamma^2;
        end
        loop_up=(dist(k)>=rho_n);
     end

end % GOTO LOOP_LABEL
