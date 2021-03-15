% Copyright 2017-2018 CNRS-UM LIRMM
% Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

function [ w ] = generate_weights( m, t, n, s )
%generate_weights Generate weights for filtering at the end of window

w = [];
for i=-m:m
w = [w Weight(i, t, m, 2, 0)];
end
end

