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


function [ res ] = savitzky_golay( x, sg )
%savitzky_golay Apply Gram Convolution weights to x
%   sg: weights
%    x: data
  res = 0;
  for i=1:length(sg)
    res = res + sg(i) * x(i);
  end
end

