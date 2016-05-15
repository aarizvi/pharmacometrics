function [ output ] = ssort( X )
%Returns vector Y in ascending order
sym Y;

%IMPLEMENTATION OF BUBBLE ALGORITHM FOR ASCENDED SORTING
n = length(X);
% n-1 passes
if size(X,1) > 1 & size(X,2) > 1 %unfortunately could not get this work for more than one dimension
    warning('ssort() only sorts 1-dimensional vectors')
else
    for j=1:1:n-1
    % comparing each adjacent number with the next and swap
        for i=1:1:n-1
            if X(i)>X(i+1); %if pair is out of order
                % swap them and remember them using 'temp'
                temp=X(i);
                X(i)=X(i+1);
                X(i+1)=temp;
            end
        end
    end
Y = X;
display(Y)    
end
end

