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

function [ gf ] = GenFact( a, b )
%GenFact Calculates the generalized factorial (a)(a-1)(a-2)...(a-b+1)
%   Detailed explanation goes here
    gf = 1;
    for j=(a-b+1):a
        gf = gf * j;
    end
end

