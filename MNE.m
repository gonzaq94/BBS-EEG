function [s]=MNE(x,A,lambda)
    
    sizeA = size(A);
    s = transpose(A)*inv(A*transpose(A)+lambda*eye(sizeA(1)))*x;