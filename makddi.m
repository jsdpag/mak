
function  [ ddi , bci ] = makddi (  r  ,  d  ,  a  )
% 
% ddi = makddi (  r  ,  d  ,  a  )
% 
% MET Analysis Kit. Implements a Disparity Discrimination Index ( DDI ).
% Takes input matrix r of neural firing rates with N trials ( rows ) and C
% spike clusters ( columns ) and N element vector d of disparities on each
% trial. Alternatively, r can be a number of rows by number of time bins
% array of firing rates from a single cluster. Returns ddi, a 1 x C DDI
% measurement for each spike cluster. Optionally returns bci, a 2 x C
% matrix of bootstrap confidence intervals around each DDI measurement
% ( 2000 bootstrap samples ). Optional input a sets the alpha value for the
% bootstrap CI ; this is 0.05 by default for 95% CI.
% 
% Take special note that firing rates in r must all be computed with the
% same length of analysis window. Analysis window length affects the SSE
% term, hence varying it renders data incomparable. Raw firing rates are
% expected, and the square root of these are first computed by the function
% before deriving DDI.
% 
% 
% Implements the equation:
%   
%   ( Rmax  -  Rmin )  /  ( Rmax  -  Rmin  +  2 * sqrt( SSE / ( N - M ) ) )
% 
% Where Rmax and Rmin are the maximum and minimum average firing rate of
% the measured tuning curve. SSE is the the square root of the residual
% variance around the mean response to each disparity. N is the number of
% trials and M is the number of measured disparities. Hence, the number of
% trials must exceed the number of measured disparities by at least 1.
% 
% NaN values in r are ignored. If there is missing data for a column of r,
% then 0 is returned rather than NaN , the assumption is that missing data
% denotes a lack of disparity discriminability.
% 
% 
% References:
%   
%   Prince SJ, Pointon AD, Cumming BG, Parker AJ. Quantitative analysis of
%     the responses of V1 neurons to horizontal disparity in dynamic
%     random-dot stereograms. J Neurophysiol. 2002 Jan;87(1):191-208.
%   
%   Uka T, DeAngelis GC. Contribution of middle temporal area to coarse
%     depth discrimination: comparison of neuronal and psychophysical
%     sensitivity. J Neurosci. 2003 Apr 15;23(8):3515-30.
%   
%   Watanabe M, Tanaka H, Uka T, Fujita I. Disparity-selective neurons in
%     area V4 of macaque monkeys. J Neurophysiol. 2002 Apr;87(4):1960-73.
% 
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% Input checking %%%
  
  % r is a numeric, real-valued matrix
  if  ~ ismatrix(  r  )  ||  ~ isnumeric(  r  )  ||  ~ isreal(  r  )  ||...
      any( isinf(  r( : )  ) )
    
    error (  'MAK:makddi:inputr'  ,  ...
      'makddi: r must be a finite, real, numeric matrix'  )
    
  % d must be a numeric, real, finite vector
  elseif  ~ isvector(  d  )  ||  ~ isnumeric(  d  )  ||  ...
      ~ isreal(  d  )  ||  ~ all( isfinite(  d  ) )
    
    error (  'MAK:makddi:inputd'  ,  ...
      'makddi: d must be a finite (no Inf or NaN), real, numeric vector'  )
    
  % Does d have same length as dimension 1 in r?
  elseif  length (  d  )  ~=  size (  r  ,  1  )
    
    error (  'MAK:makddi:dlength'  ,  ...
      'makddi: length of d must match length of dim 1 in r'  )
    
  % Alpha given
  elseif  2  <  nargin
    
    % Check it is scalar, numeric, real, finite, and between 0 and 1
    if  ~ isscalar (  a  )  ||  ~ isnumeric (  a  )  ||  ...
        ~ isreal (  a  )  ||  ~ isfinite (  a  )  ||  a < 0  ||  1 < a
      
      error (  'MAK:makddi:inputa'  ,  ...
      'makddi: a must be a finite, real, numeric scalar between 0 and 1'  )
      
    end
    
  % No alpha given
  else
    
    % Set default
    a = 0.05 ;
    
  end % check input
  
  
  %%% Prep %%%
  
  % Guarantee that d is a column vector
  if  ~ iscolumn (  d  )  ,  d = d( : ) ;  end
  
  % Take the square root of firing rates
  r = sqrt (  r  ) ;
  
  % Get grouping information , grouping trials by measured disparity.
  % Because there is one grouping term , the number of columns in G will
  % equal the number of unique disparities.
  G = makfun (  d  ) ;
  
  % Check that the number of trials exceeds the number of disparity levels
  if  size (  G  ,  1  )  <=  size (  G  ,  2  )
    
    error (  'MAK:makddi:ntrials'  ,  [ 'makddi: the number of ' , ...
      'trials must exceed the number of disparity levels' ]  )
    
  end % check num trials
  
  
  %%% Disparity discrimination index %%%
  
  % Empirical value
  ddi = fddi (  r  ,  G  ) ;
  
  % No confidence intervals requested , skip bootstrapping
  if  nargout  <  2  ,  return  ,  end
  
  % Bootstrap confidence intervals
  parfor  p = 1 : size (  r  ,  2  )
    
    bci( : , p ) = bootci (  2e3  ,  { @fddi , r( : , p ) , G }  ,  ...
      'alpha'  ,  a  ) ;
    
  end
  
  
end % makddi


%%% Subroutines %%%

function  ddi = fddi (  r  ,  G  )
  
  % Total number of trials and disparities
  [ ~ , M ] = size (  G  ) ;
  
  % Calculate disparity tuning curve , return disparity by cluster matrix
  c = makfun (  @( x ) mean( x , 1 , 'omitnan' )'  ,  r  ,  G  ,  M  )' ;
  
  % Guarantee that disparities indexed along rows if r has one column
  if  isvector (  r  )  ,  c = c' ;  end
  
  % Get maximum and minimum responses
  Rmax = max (  c  ,  []  ,  1  ) ;
  Rmin = min (  c  ,  []  ,  1  ) ;
  
  % Take the difference
  Rdif = Rmax  -  Rmin ;
  
  % Initialise sum of squared errors to zero
  SSE = zeros ( size(  Rmax  ) ) ;
  
  % Disparity levels
  for  i = 1 : M
    
    % Set of trials at this disparity level
    j = G( : , i ) ;
    
    % Accumulate SSE at this disparity level
    SSE( : ) = SSE  +  sum (  ...
      bsxfun( @minus , r( j , : ) , c( i , : ) )  .^  2  ,  ...
        1  ,  'omitnan'  ) ;
    
  end % disparities
  
  % Count number of non-nan values for each cluster
  N = sum (  ~isnan(  r  )  ,  1  ) ;
  
  % Disparity discrimination index denominator
  ddi = Rdif  +  2 * sqrt(  SSE  ./  ( N - M )  ) ;
  
  % Disparity discrimination index
  ddi( : ) = Rdif  ./  ddi ;
  
  % Set any columns where N <= M to NaN
  ddi(  N  <=  M  ) = NaN ;
  
  % Replace NaN values with zero
  ddi( isnan(  ddi  ) ) = 0 ;
  
end % fddi

