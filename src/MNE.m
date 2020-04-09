function [s]=MNE(x,A,lambda)
    
    s = transpose(A)*inv(A*transpose(A)+lambda*eye(size(A,1)))*x;