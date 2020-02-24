function [ s ] = SISSY(x,A,T,lambda,alpha,Niterations)
%SISSY Summary of this function goes here
%   Detailed explanation goes here
    rho = 1;
    
    P=sparse(rho*(T.'*T+speye(size(T,2))));
    APi=A/P;
    L=chol(eye(size(A,1))+APi*A.','lower');
    s1=A.'*x;

    z = zeros(size(T,1),size(x,2));
    u = zeros(size(T,1),size(x,2));
    y = zeros(size(T,2),size(x,2));
    v = zeros(size(T,2),size(x,2));
    s = zeros(size(T,2),size(x,2));
    
    for it=1:Niterations
        
        b=s1+rho*(T.'*(z+u/rho)+y+v/rho);
        s=P\b-APi'*(L'\(L\(APi*b)));
        
        z = prox(T*s-u/rho,lambda/rho);
        
        y = prox(s-v/rho,lambda*alpha/rho);
        
        u = u + rho*(z-T*s);
        
        v = v + rho*(y-s);
        
    end
    
end

