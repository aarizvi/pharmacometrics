function [ output ] = mmax( X )
%Maximum of n-dimensional vector without using max() function
if size(X,1) > 1
    ans = sort(X, 'descend');
    frow = ans(1,:);
    frow = sort(frow, 'descend');
    max = frow(1);
    max
else
    ans = sort(X, 'descend');
    max = ans(1);
    max
end

