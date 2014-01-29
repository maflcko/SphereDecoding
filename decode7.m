function u_closest = decode7( y, H )
%DECODE7( y, H ) Lattice Decode
%   u_closest = DECODE7( y, H ) returns the indexes of the closest lattice
%   vector for a given inverse lower triangular generator matrix H and any
%   row vector y.
%   (Algorithm 7 by A. Ghasemmehdi, E. Agrell)

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
k=n+1;
d=zeros(1,n)+n;
dist(n+1)=0;
E(n,:)=y*H;

gamma(n)=0; %init needed due to matlab 'restrictions'

while(true) %LOOP_LABEL
    
    loop_down=true;
    while(loop_down)
        if(~(k==1))
            k=k-1;
            for a=d(k):-1:(k+1)
                E(a-1,k)=E(a,k)-gamma(a)*H(a,k);
            end
            u(k)=round(E(k,k));
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
    loop_up=true;
    while (loop_up)
        if (k==n)
            return;
        else
            k=k+1;
            u(k)=u(k)+step(k);
            step(k)=-step(k)-sgn(step(k));
            gamma(k)=(E(k,k)-u(k))/H(k,k);
            dist(k)=dist(k+1)+gamma(k)^2;
        end
        loop_up=(dist(k)>=rho_n);
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
