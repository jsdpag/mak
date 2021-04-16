
function  rp = makrpcorr (  varargin  )
% 
% rp = makrpcorr (  ...  )
% 
% MET Analysis Kit. A wrapper function for Matlab's corr ( ). Packs the RHO
% and PVAL outputs into a single N by 2 matrix, with RHO and PVAL in
% columns 1 and 2. Use reshape on one of the columns to restore the
% original RHO or PVAL as returned by corr.
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  % Run corr
  [ r , p ] = corr (  varargin{ : }  ) ;
  
  % Pack output
  rp = [  r( : )  ,  p( : )  ] ;
  
end % makrpcorr

