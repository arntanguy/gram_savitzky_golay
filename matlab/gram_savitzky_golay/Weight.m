% Copyright 2017-2018 CNRS-UM LIRMM
% Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

function [ w ] = Weight( i, t, m, n, s )
%Weight Calculates the weight of the ith data point for the t'th
%Least-Square point of the s'th derivative, over 2m+1 points, order n

w = 0;
for k=0:n
    w = w + (2*k+1) * (GenFact(2*m,k)/GenFact(2*m+k+1,k+1))*GramPoly(i,m,k,0)*GramPoly(t,m,k,s);
end
end

