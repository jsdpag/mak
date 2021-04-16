
function  E = makenergymat (  n  ,  ca  ,  c  ,  d0  )
% 
% E = makenergymat (  n  ,  ca  ,  c  ,  d0  )
%
% MET Analysis Kit, pre-processing. After initial clustering, the next
% step in spike sorting (Fee et al. 1996) is to compute the interface-
% energy matrix, which tabulates the interface energy between every pair of
% spike clusters. The un-normalised energy matrix is returned, as this can
% be easily updated when clusters are later merged. Hence, it is necessary
% to compute connection strengths from these raw energy values.
% 
% 
% Input
% 
%   n - Vector of the number of spikes assigned to each cluster. n( i ) is
%     the number of spikes in the ith cluster.
% 
%   ca - Vector of cluster assignments for each spike. ca( i ) returns the
%     cluster assignment of the ith spike.
%   
%   c - S x N matrix of spike waveform components. Spikes are indexed
%     across columns, and components over rows. Hence c( : , i ) is the set
%     of components for the ith spike waveform.
%   
%   d0 - Scalar value , the scaling term returned by makspkclust.
% 
% 
% Output
% 
%   E - Nc x Nc matrix , where Nc is the number of initial spike clusters.
%     Values are tabulated in the upper-triangular section of the matrix
%     and along the diagonal. Thus E( i , j ), where i <= j, is the raw
%     interface energy between spike clusters i and j.
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
  
  
  %%% Preparation %%%
  
  % The number of clusters
  cnum = numel (  n  ) ;
  
  % Group spikes by cluster number , for distribution of data in parfor
  % loop
  C = arrayfun (  @( i ) c(  :  ,  ca  ==  i  )'  ,  ...
    1 : cnum  ,  'UniformOutput'  ,  false  ) ;
  
  % Make vectors of cluster data in the order that they're needed to
  % compute energies
  npairs = ( cnum ^ 2 - cnum ) / 2  +  cnum ;
  C1 = cell (  npairs  ,  1  ) ;
  C2 = cell (  npairs  ,  1  ) ;
  k = 0 ;
  
  % In this order, we fill in each row of the upper-triangular portion of
  % the energy matrix
  for  i = 1 : cnum
    for  j = i : cnum
      k = k + 1 ;
      C1{ k } = C{ i } ;
      C2{ k } = C{ j } ;
    end
  end
  
  
  %%% Raw interface Energy %%%
  
  parfor  k = 1 : npairs
    
    % Pairwise distances between all spikes in cluster i with those in
    % cluster j
    d = pdist2 (  C1{ k }  ,  C2{ k }  ) ;
    
    % Compute interface energy
    e( k ) = sum ( exp(  - double( d( : ) )  /  d0  ) ) ;
    
  end
  
  % Reshape from a vector into a square matrix , first allocate a matrix
  E = zeros (  cnum  ) ;
  
  % Then get logical index that fills lower triangular portion. Since
  % Matlab indexing is columns first, we will transpose the result to get
  % the final energy matrix.
  i = tril ( true(  cnum  ) ) ;
  
  % Assign values
  E( i ) = e ;
  
  % And transpose
  E = E' ;
  
  
  %%% Correction %%%

  % Index vector of diagonal elements
  i = find ( eye(  cnum  ) ) ;

  % Correction of terms , following UltraMegaSort2000's recipe in ss_energy
  E( i ) = ( E( i )  -  double( n' ) )  /  2 ;
  
  
end % makenergymat

