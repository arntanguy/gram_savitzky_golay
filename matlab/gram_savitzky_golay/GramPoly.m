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


function [ gp ] = GramPoly( i, m, k, s )
%GRAMPOLY Calculates the Gram Polynomial (s=0) or its sth derivative 
% evaluated at i, order k, over 2m+1 points

gp = 0;
if(k>0)
    gp = (4*k-2)/(k*(2*m-k+1)) * (i* GramPoly(i,m,k-1,s) + ...
    s * GramPoly(i,m,k-1,s-1)) - ((k-1)*(2*m+k))/(k*(2*m-k+1))*GramPoly(i,m,k-2,s);
else
    if (k==0) && (s==0)
        gp = 1;
    else
        gp = 0;
    end
end

end

