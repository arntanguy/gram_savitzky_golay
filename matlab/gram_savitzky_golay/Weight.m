% Copyright 2017-2018 CNRS-UM LIRMM
% Copyright 2017-2018 Arnaud TANGUY <arnaud.tanguy@lirmm.fr>
%
% This file is part of robcalib.
%
% robcalib is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% robcalib is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with robcalib.  If not, see <http://www.gnu.org/licenses/>.

function [ w ] = Weight( i, t, m, n, s )
%Weight Calculates the weight of the ith data point for the t'th
%Least-Square point of the s'th derivative, over 2m+1 points, order n

w = 0;
for k=0:n
    w = w + (2*k+1) * (GenFact(2*m,k)/GenFact(2*m+k+1,k+1))*GramPoly(i,m,k,0)*GramPoly(t,m,k,s);
end
end

