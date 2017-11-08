function [ gf ] = GenFact( a, b )
%GenFact Calculates the generalized factorial (a)(a-1)(a-2)...(a-b+1)
%   Detailed explanation goes here
    gf = 1;
    for j=(a-b+1):a
        gf = gf * j;
    end
end

