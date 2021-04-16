
function  p = makpca (  pv  ,  w  )
% 
% p = makpca (  pv  ,  w  )
% 
% MAK Analysis Kit, pre-processing. Performs principal component analysys
% via singular value decomposition the a set of waveforms in w. The first N
% components of each waveform are returned in p ; N is the number of
% components required to capture pv percent of the variance.
% 
% Input
% 
%   pv - The minimum percentage of variance to be captured by the first N
%     components returned from each waveform.
% 
%   w - S x M matrix of M waveforms consisting of S samples each.
% 
% Output
% 
%   p - N x M matrix of M waveforms consisting of the first N components
% 
% 
% See:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Disable warning message:
  % Columns of X are linearly dependent to within machine precision.
  % Using only the first 26 components to compute TSQUARED.
  wstate = warning (  'off'  ,  'stats:pca:ColRankDefX'  ) ;
  
  % Perform principal component analysis on all aligned waveforms
  [ coef , ~ , ~ , ~ , explained ] = pca (  w'  ) ;
  
  % Restore previous warning state
  warning (  wstate  ) ;
  
  % Number of components to keep , in descending order of component
  % variance
  n = find (  pv  <  cumsum( explained )  ,  1  ,  'first'  ) ;

  % Throw away unwanted coefficients
  coef( : , n + 1 : end ) = [] ;

  % Compute component of each aligned waveform
  p = coef'  *  w ;
  
end % makpca

