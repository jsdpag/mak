
function  [ n , c , E , m ] = makcmerge (  n  ,  c  ,  E  ,  cut  )
% 
% [ n , c , E , m ] = makcmerge (  n  ,  c  ,  E  ,  cut  )
% 
% MET Analysis Kit, pre-processing. Initial spike clusters are merged
% together on the basis of their connection strength. Merging is done
% following these steps:
% 
%   1) Find pair of clusters with the highest connection strength
%   2) Merge clusters and compute new raw interface energies
%   3) Repeat 1 and 2 until a cutoff connection strength is reached
% 
% 
% Input
% 
%   n - A row vector containing the number of spikes per cluster
%   
%   c - A row vector containing the cluster assignment of each spike
%   
%   E - Raw interface-energy matrix
%   
%   cut - Connection strength cutoff value
% 
% 
% Output
% 
%   n - Final number of spikes per cluster
%   
%   c - Final cluster assignment of all spikes
%   
%   E - Final raw interface energy matrix
%   
%   m - 2 x M list of cluster mergers. The top row holds the lower cluster
%     number of each merger, and the bottom row holds the higher. Mergers
%     are ordered as they happened from left-to-right.
% 
% 
% References:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% UltraMegaSort2000, https://neurophysics.ucsd.edu/software.php
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Number of clusters
  cnum = numel (  n  ) ;
  
  % There is only one cluster , return empty merger list
  if  cnum  ==  1
    m = zeros (  2  ,  0  ,  'uint8'  ) ;
    return
  end
  
  % Table of untested pairs. Clusters cannot be compared with themselves.
  U = triu (  true( cnum )  ,  1  ) ;
  
  % Connection strength mask
  M = ~ U ;
  
  % Allocate list of cluster mergers
  mi = 0 ;
  m = zeros (  2  ,  cnum - 1  ,  'uint8'  ) ;

  % Aggregation loop runs so long as there are untested pairs
  while  any (  U( : )  )

    % Compute new connection strengths
    J = makconnstrength (  E  ,  n  ) ;

    % Mask tested connection strengths
    J( M ) = 0 ;

    % Find cluster pair with the maximum connection strength
    [ cmax , pmax ] = max (  J( : )  ) ;

    % Connection strengths have fallen to below the cutoff point , end now
    if  cmax  <  cut  ,  break  ,  end
    
    % Get the number of each cluster in the pair. c1 and c2 are necessarily
    % the low and high cluster numbers because of the upper-triangular
    % shape of E.
    [ c1 , c2 ] = ind2sub (  [ cnum , cnum ]  ,  pmax  ) ;

    % Mark that this pair is no longer un-tested
    U( c1 , c2 ) = 0 ;
    M( c1 , c2 ) = 1 ;
    
    % Record the merger
    mi = mi  +  1 ;
    m( : , mi ) = [  c1  ;  c2  ] ;
    
    % Re-assign spikes in the high-number cluster to the low-number cluster
    c( c  ==  c2 ) = c1 ;

    % Compute new intra-cluster energy
    E( c1 , c1 ) = E( c1 , c1 )  +  E( c2 , c2 )  +  E( c1 , c2 ) ;

    % Get linear indices of for inter-cluster energies between the low
    % cluster and all others -- with the exception of the high cluster.
    % Takes an 'L' shape out of the upper triangle of the raw energy table.
    [ i1 , i2 , j ] = makcind (  c1  ,  c2  ,  cnum  ) ;

    % Compute new inter-cluster energies
    E( i1 ) = E( i1 )  +  E( i2 ) ;

    % Mark that the high-number cluster no longer needs testing with any
    % other pair. Also zero energy matrix for inter-cluster energies and
    % set intra-cluster energy to 1, according to UltraMegaSort2000
    % convention.
    U( i2 ) = 0 ;
    M( i2 ) = 1 ;
    E( i2 ) = 0 ;
    E( c1 , c2 ) = 0 ;
    E( c2 , c2 ) = 1 ;

    % Update number of waveforms in new cluster
    n( c1 ) = n( c1 )  +  n( c2 ) ;

    % Old cluster is empty
    n( c2 ) = 0 ;

    % We can now make sure that the new cluster is compared with any
    % non-empty cluster
    U( i1(  0  <  n( j )  ) ) = 1 ;

  end % agg loop
  
  % Remove empty tail from merge list
  m( : , mi + 1 : end ) = [] ;
  
  % Convert interface energy matrix to sparse matrix
  E = sparse (  E  ) ;
  
end % makcmerge

