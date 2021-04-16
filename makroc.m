
function  [ auc , varargout ] = makroc (  x  ,  p  ,  a  )
% 
% [ auc , y , ci , t , f , th ] = makroc (  x  ,  p  ,  a  )
% 
% MET Analysis Kit. Computes the Receiver Operating Characteristic curves
% for data in x and optionally returns the area under each curve ( auc ),
% the Youden J statistic ( y ) , the bootstrap confidence interval around
% auc ( ci ), the true positive rate ( t ), the false positive rate ( f ),
% and the threshold level ( th ).
% 
% 
% x must be a real-valued ( i.e. no imaginary component ) numeric matrix of
% N rows with any size in dimensions 2 and higher. It is considered to be N
% samples of multivariate random variable X. Input p is a logical vector of
% N elements that is true ( 1 ) for every sample ( row ) of x that is a
% true positive, and false ( 0 ) for every sample of x that is a false
% positive. A separate ROC curve is computed for each column of x . That
% is, a ROC curve is computed for each univariate marginal of random
% variable X. No Inf or NaN values are allowed.
% 
% In analysing neural responses, x could be the firing rate of M
% simultaneously recorded units over N trials. Or it could be the firing
% rate of one unit for T time bins over N trials. It could even be the
% firing rate on N trials over T time bins for M units , arranged in a 3D
% matrix. If p identifies trials where one choice was made versus another
% then auc can be interpreted as the choice probability for each unit /
% time bin.
% 
% a is an optional input that sets the alpha value, or level of
% significance at which bootstrap intervals in ci are computed. This
% defaults to 0.05. If an empty matrix [] is provided then bootstrap
% intervals are not computed and ci will return an empty matrix. This may
% be desirable when outputs t, f, and th are required but bootstrap
% intervals are not.
%
% 
% auc contains areas under the ROC curve for each column in x. It's size
% matches that of x starting with dimension 2.
% 
% Youden's J statistic is the threshold level that provides the maximum
% discrimination between true and false positives. In this case, it is
% taken to be the threshold that maximises the absolute distance from the
% ROC curve to the chance line. Specifically, it gives the threshold of the
% true positive rate that is the furthest above or below the chance line.
% y matches auc in size.
% 
% ci gives the bootstrap confidence intervals around each value in auc. It
% has 2 rows and otherwise matches x in size. Row 1 contains the lower
% confidence interval and row 2 the upper.
% 
% t, f, and th give the true ( t ) and false ( f ) positive rates for each
% observed threshold value ( th ) in x. Each matches x in size , plus one
% extra row for threshold value - Inf. Hence one may plot the true or false
% positive rate as a function of threshold by plotting th on the x-axis
% against values of t or f in the y-axis. Alternatively, the ROC curves can
% be plotted when false positive rates are on the x-axis and true positive
% rates are on the y-axis. t , f , and th values are all set to NaN for
% repeats of any given threshold value ; only the final instance ( row ) of
% each unique threshold value contains non-NaN data.
% 
% 
% NOTE : Requires Matlab's Parallel Processing Toolbox to compute bootstrap
%   intervals
% 
% 
% References:
%  
%   Britten KH1, Newsome WT, Shadlen MN, Celebrini S, Movshon JA. A
%     relationship between behavioral choice and the visual responses of
%     neurons in macaque MT. Vis Neurosci. 1996 Jan-Feb;13(1):87-100.
%   
%   Green DM, Swets JA. (1966). Signal detection theory and psychophysics.
%     New York, NY: John Wiley and Sons Inc. ISBN 0-471-32420-5. 
%   
%   Youden WJ. Index for rating diagnostic tests. Cancer. 1950
%     Jan;3(1):32-5.
% 
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Check maximum and minimum number of input/output arguments
   narginchk (  2  ,  3  )
  nargoutchk (  0  ,  6  )
  
  % Size of x
  s = size (  x  ) ;
  
  % x must be a real-valued numeric matrix
  if  ~ isnumeric (  x  )  ||  ~ isreal (  x  )  ||  ...
      ~ all ( isfinite(  x( : )  ) )
    
    error (  'MAK:makroc:x'  ,  ...
      'makroc: x must be a real, numeric matrix of finite values'  )
    
  % p must be a logical vector with as many elements as there are rows in x
  elseif  ~ isvector (  p  )  ||  ~ islogical (  p  )  ||  ...
      numel (  p  )  ~=  s( 1 )
    
    error (  'MAK:makroc:p'  ,  [  'makroc: p must be a logical ' , ...
      'vector with length equal to the size of dim 1 in x' ]  )
    
  % a given
  elseif  2  <  nargin
    
    % Allowed to be empty
    if  isempty (  a  )
      
      % No action needed
      
    % Check that a is scalar numeric value between 0 and 1
    elseif  ~ isscalar (  a  )  ||  ~ isnumeric (  a  )  ||  ...
        ~ isreal (  a  )  ||  a  <  0  ||  1  <  a
      
      error (  'MAK:makroc:a'  ,  [  'makroc: a must be a scalar, ' , ...
      'real, numeric value between 0 and 1' ]  )
      
    end % Check a
    
  % a not given , set to default of 0.05
  else
    
    a = 0.05 ;
    
  end % check input
  
  % But x is empty to return empties now
  if  isempty (  x  )
    auc = [] ;
    for  i = 2 : nargout  ,  varargout{ i } = [] ;  end %#ok
    return
  end
  
  
  %%% Prepare data %%%
  
  % Number of rows and columns
  N = s( 1 ) ;
  M = prod (  s( 2 : end )  ) ;
  
  % Unravel x into a 2D matrix
  x = reshape (  x  ,  N  ,  M  ) ;
  
  % Sort values ascending in each column of x
  [ x , i ] = sort (  x  ,  1  ) ;
  
  % Replicate p across all columns of x
  p = repmat (  p( : )  ,  1  ,  M  ) ;
  
  % Sort each column of p to match the order in each column of x
  p = p(  i  ) ;
  
  % And have p match the type of x if x is not double
  if  ~ isa (  x  ,  'double'  )
    
    p = cast (  p  ,  'like'  ,  x  ) ;
    
  end % x not double
  
  
  %%% Empirical results %%%
  
  [ auc , t , f ] = roc (  x  ,  p  ,  []  ) ;
  
  % Match auc to original shape of x
  auc = reshape (  auc  ,  [ s( 2 : end ) , 1 ]  ) ;
  
  % Return Youden's J statistic
  if  1  <  nargout
    
    % Notice that the fraction of false positives that provide the ROC
    % curve x-axis values is also the height of the chance line on the y-
    % axis , because we're dealing with a right-angle triangle where the
    % hypotenuse lines up with the chance line. Therefore, locate the index
    % in each column where the difference between true and false positive
    % fraction is maximised.
    [ ~ , i ] = max (  abs( t  -  f )  ,  []  ,  1  ) ;
    
    % x has one less row than t or f , so subtract 1 to bring i into
    % register with x
    i = i  -  1 ;
    
    % However , we may have some zeros in there now. This could happen if
    % the true and false-positive distributions were identical , then t
    % minus f will be zero at all thresholds and the index returned by max
    % will be 1. Remember where the zeros were and reset them to 1 to avoid
    % an indexing error.
    i_zero = i  <  1 ;
    i( i_zero ) = 1 ;
    
    % Get subscripts for columns
    I = 1 : M ;
    
    % Now we can convert multidimensional subscripts to linear indices
    i = sub2ind (  s  ,  i  ,  I  ) ;
    
    % Keep column subscripts for cases where t - f = 0
%     I = I( i_zero ) ;
%     
%       % Guarantee that I has one row
%       I = reshape (  I  ,  1  ,  numel( I )  ) ;
%     
%     % Since these cases were reset to a row index of one , we can use a
%     % vector of 1's for dim 1 sub-scripts
%     i_zero = ones (  1  ,  sum(  i_zero( : )  )  ) ;
%     
%     % And get linear indices to those cases
%     i_zero = sub2ind (  s  ,  i_zero  ,  I  ) ;
    
    % Get threshold for these points , and match the size of dims 2 and
    % higher in x from dim 1 in the output. That is , we squeeze the size
    % of x down one dimension for Youden J statistics as there are no
    % trials along rows any more.
    varargout{ 1 } = x(  i  ) ;
    
    % But remember the cases where t - f = 0? The row index returned by max
    % was 1 in register with t and f , which is the proportion of true- and
    % false-positives greater than the value of - Inf. So return - Inf.
    varargout{ 1 }( i_zero ) = - Inf ;
    
    % Reshape output to match original shape of x
    varargout{ 1 } = reshape (  varargout{ 1 }  ,  [ s( 2 : end ) , 1 ]  );
    
  end % Youden's J
  
  % Return true positive rate
  if  3  <  nargout
    
    % Match original shape of x
    t = reshape (  t  ,  [ N + 1 , s( 2 : end ) ]  ) ;
    
    % Return
    varargout{ 3 } = t ;
  
  end % TPR
  
  % Return false positive rate
  if  4  <  nargout
    
    % Match original shape of x
    f = reshape (  f  ,  [ N + 1 , s( 2 : end ) ]  ) ;
    
    % Return
    varargout{ 4 } = f ;
  
  end % FPR
  
  % Return threshold values
  if  5  <  nargout
    
    % - Inf is always the first threshold value
    varargout{ 5 } = [  -inf(  1  ,  M  )  ;  x  ] ;
    
    % Match original shape of x
    varargout{ 5 } = ...
      reshape (  varargout{ 5 }  ,  [ N + 1 , s( 2 : end ) ]  ) ;
    
    % Set any value to NaN that is not the final instance of each unique
    % threshold
    varargout{ 5 }(  isnan( t )  ) = NaN ;
  
  end % ret thresh
  
  
  %%% Bootstrap confidence intervals %%%
  
  % ci output requested
  if  2  <  nargout
    
    % But no boot ci is desired
    if  isempty (  a  )
      
      % Return empty , then
      varargout{ 2 } = [] ;
      
    % Boot ci desired
    else
      
      % Halve alpha value so that we can get the lower and upper interval
      % embracing 100 * ( 1 - alpha ) percent of the area under the
      % bootstrap sample distribution
      a = a  /  2 ;
      
      % Generate boot samples
      parfor  i = 1 : 2e3
        
        % Randomly sample trials with replacement
        j = ceil (  N  *  rand (  N  ,  1  )  ) ;
        
        % However , we must sort trial numbers so that we don't need to
        % re-sort x or p. This saves time for large data sets.
        j = sort (  j  ) ;
        
        % Compute area under curve , beware that Matlab editor complains
        % about how x and p are now broadcast variables , but how else do
        % we do we return a sorted sub-set?
        aboot( i , : ) = roc (  x( j , : )  ,  p( j , : )  ) ; %#ok
        
      end % boot samples
      
      % Allocate return confidence intervals
      ci = zeros (  [ 2 , s( 2 : end ) ]  ,  'like'  ,  x  ) ;
      
      % Get percentiles for confidence interval levels
      ci( : ) = prctile (  aboot  ,  100 * [ a , 1 - a ]  ,  1  ) ;
      
      % Return confidence intervals
      varargout{ 2 } = ci ;
    
    end % compute boot ci?
    
  end % boot ci
  
  
end % makroc


%%%  Subroutines  %%%

% Calculate ROC curves from data in x with classifications in p. Returns
% the true and false positive rates. t and f are the size of x plus 1 extra
% row. Their rows are in register with all possible ROC threshold values ,
% including repeats. NaN is set for all rows that do not contain the final
% instance of each unique threshold value. Provide dummy 3rd input arg , an
% empty matrix say , to signal that an optimisation for data with no repeat
% values can be tried.
function [ auc , t , f ] = roc (  x  ,  p  ,  varargin  )
  

  %-- Preparation --%
  
  % Number of rows and columns in the data
  [ N , M ] = size (  x  ) ;
  
  % Allocate logical array. These will be in register with values of t and
  % f that store the complete true and false positive rates counted over
  % repetitions of each threshold value. First row is the unique instance
  % of threshold value - Inf.
  uni = true (  1 + N  ,  M  ) ;
  
  % Find the last instance of each unique threshold value in x
  uni( 2 : end - 1 , : ) = x( 1 : end - 1 , : )  ~=  x( 2 : end , : )  ;
  
  % Find cases where there are no repeat threshold values
  if  2  <  nargin  ,  norep = all (  uni  ,  1  ) ;  end
  
  
  %-- True and false positive rates --%
  
  % Number of true positives and false positives
  Nt = sum (  p  ,  1  ) ;
  Nf = N  -  Nt ;
  
  % Count the number of true and false positives for each observation ,
  % normalise to fractions
  t = bsxfun (  @rdivide  ,  cumsum(   p  ,  1  )  ,  Nt  ) ;
  f = bsxfun (  @rdivide  ,  cumsum( ~ p  ,  1  )  ,  Nf  ) ;
  
  % But of course, we now have the empirical cumulative distribution
  % functions that say how many values are less than or equal to each
  % possible observation. Subtract from 1 to get the fraction of values
  % above each possible threshold. Append zeros to the top of each CDF ,
  % first.
  t = 1  -  [  zeros(  1  ,  M  )  ;  t  ] ;
  f = 1  -  [  zeros(  1  ,  M  )  ;  f  ] ;
  
  
  %-- Area under ROC --%
  
  % Allocate return array
  auc = zeros (  1  ,  M  ,  'like'  ,  x  ) ;
  
  % Find cases where there are repeat threshold values
  if  2  <  nargin
    
    % Signal arg given , get linear index of cases with repeats. Guarantee
    % that a row vector of indices or empty matrix without non-zero
    % dimensions is returned.
    i = find (  ~ norep  ) ;
    if  ~ isrow (  i  )  ,  i = i( : )' ;  end
    if  isempty (  i  )  ,  i = [] ;  end
    
  else
    
    % No signal arg given , all cases are assumed to have repeats , such as
    % when bootstrapping is done
    i = 1 : M ;
    
  end % find cases with repeats
  
  
  % Columns , we will use a combination of subscript indexing to access
  % rows and linear indexing to access columns spanning dimensions 2 and
  % higher.
  for  i = i
    
    % Get subscript indices for rows containing the last instance of each
    % unique threshold value
    j = find (  uni( : , i )  ) ;
    
    % j and k line up each successive pair of true-/false-positive rates
    k = j( 2 : end ) ;
    j = j( 1 : end - 1 ) ;
    
    % Perform a trapezoidal integration of the curves , see dot trapz tips.
    % This is not yet normalised with division by two or corrected with a
    % sign change ... see below.
    auc( i )  =  sum (  ( f( k , i )  -  f( j , i ) )  .*  ...
                        ( t( k , i )  +  t( j , i ) )  ) ;
    
  end % columns
  
  % Signal arg given , attempt optimisation for data with no repeat values
  if  2  <  nargin
    
    % Perform integration for all such cases
    auc( norep ) = ...
    sum (  ( f( 2 : end , norep )  -  f( 1 : end - 1 , norep ) )  .*  ...
           ( t( 2 : end , norep )  +  t( 1 : end - 1 , norep ) )  ,  1  ) ;
                  
  end % optimisation
  
  % The unary negation accounts for the fact that f and t start at 1 and
  % work towards 0 , resulting in a negative area. If we reversed the
  % order of f and t then we would get the same absolute value for area,
  % but positive. Divide by 2 to complete trapezoidal integration.
  auc = - auc  /  2 ;
  
  
  %-- Set t and f elements to NaN --%
  
  % t and f not requested , quit here
  if  nargout  <  2  ,  return  ,  end
  
  % Find all rows that are not in register with the final instance of each
  % unique threshold value
  uni = ~ uni ;
  
  % And set corresponding elements of t and f to NaN
  t( uni ) = NaN ;
  f( uni ) = NaN ;
  
  
end % roc

