function [ x ] = InverseGamma( a, b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x = norm(randn(round(a*2),1))^2/(2*b);

end

