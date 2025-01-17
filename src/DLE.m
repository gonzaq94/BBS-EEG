function [ DLE ] = DLE(idx,idx_est,r_grid)
%DLE Summary of this function goes here
%   Detailed explanation goes here

    %if (length(idx_est)==0)
     %   DLE = 0;
    %else

        r_idx = repmat(r_grid(idx,:),[1,1,size(r_grid(idx_est),1)]);
        r_idx_est = repmat(r_grid(idx_est,:),[1,1,size(r_grid(idx),1)]);

        s1 = mean(min(sqrt(sum((r_idx-permute(r_idx_est,[3 2 1])).^2,2)),[],1));
        s2 = mean(min(sqrt(sum((r_idx-permute(r_idx_est,[3 2 1])).^2,2)),[],3));

        DLE = (s1+s2)/2;

end

