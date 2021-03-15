% Copyright 2017-2018 CNRS-UM LIRMM
% Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

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

