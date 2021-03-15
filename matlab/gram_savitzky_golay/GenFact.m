% Copyright 2017-2018 CNRS-UM LIRMM
% Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

function [ gf ] = GenFact( a, b )
%GenFact Calculates the generalized factorial (a)(a-1)(a-2)...(a-b+1)
%   Detailed explanation goes here
    gf = 1;
    for j=(a-b+1):a
        gf = gf * j;
    end
end

