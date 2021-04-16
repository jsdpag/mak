
function  c = makcutoff (  p  ,  E  ,  n  )
% 
% c = makcutoff (  p  ,  E  ,  n  )
% 
% MET Analysis Kit, pre-processing. Estimates connection strength cutoff
% point for cluster merging by returning the upper bca bootstrap
% confidence interval of a given percentile of the connection strength
% between different clusters.
% 
% The raw interface-energy matrix in E and the number of spikes per cluster
% in n are used to compute connection strength. p is a parameter struct
% with fields nboot, alpha, and ptile. nboot bootstrap samples are taken
% and the ptile percentile is found for each, producing a distribution of
% nboot percentile values. From this, the ( 1 - alpha ) * 100% confidence
% intervals are found. If there are only two or fewer clusters then zero is
% returned.
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Number of clusters
  cnum = numel (  n  ) ;
  
  % One cluster , return
  if  cnum  <=  2
    c = 0 ;
    return
  end
  
  % Logical indices of upper triangular portion of E or J
  i = triu (  true( cnum )  ,  1  ) ;
  
  % Connection strength matrix
  J = makconnstrength (  E  ,  n  ) ;
  
  % Bootstrap confidence intervals
  ci = bootci (  p.nboot  ,  ...
    { @( x ) prctile( x , p.ptile ) , J( i ) }  ,  'alpha'  ,  p.alpha  ) ;
  
  % Return lower confidence interval
  c = max (  ci  ) ;
  
end % makcutoff

