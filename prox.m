function [ x ] = prox( y, beta )
%prox Summary of this function goes here
%   Detailed explanation goes 

    x = zeros(size(y));
    
    for i=1:length(y)
       
        if y(i) > beta
            x(i) = y(i)-beta;
        else
            if y(i) < -beta
                x(i) = y(i)+beta;
            else
                x(i) = 0;
            end
        end
    end


end

