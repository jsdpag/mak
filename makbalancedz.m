
function  [ zb , mbal , sbal ] = makbalancedz (  x  ,  g  ,  b  ,  ...
  varargin  )
% 
% [ zb , mbal , sbal ] = makbalancedz (  x  ,  g  ,  b  )
% ... = makbalancedz (  ...  ,  outtype  )
% ... = makbalancedz (  ...  ,  nanflag  )
% 
% MET Analysis Kit. Compute the balanced z-score of data in x for each
% stimulus condition defined by grouping factor g and balancing factor b.
% 
% 
% Input
% 
%   x - Column vector of N elements or multi-dimensional matrix of N rows.
%     Stores data to have balanced z-scoring. Samples are indexed across
%     rows. Columns and so forth may index any thing else , spike clusters
%     or time bins for example.
%   
%   g - N-row matrix with 1 or more columns. Each column defines a grouping
%     factor and its value for each trial. Each unique combination of
%     grouping factor values defines a sub-set of trials over which data in
%     x will undergo balanced z-scoring. For instance, g may contain the
%     baseline stimulus value, or the baseline value in column 1 and the
%     popout stimulus value in column 2.
%   
%   b - N-element logical vector. This is the balancing grouping factor.
%     That is, it identifies the two sub-sets of trials that will be used
%     to perform balanced z-scoring in each sub-set defined by g. b will be
%     true ( 1 ) for all trials of one type and false ( 0 ) for all trials
%     of another. In computing choice probability, b might identify all
%     trials in which one particular choice was made.
%   
%   outtype , nanflag - Optional string inputs instruct mean and std what
%     numeric type to return and whether to ignore NaN values. The dim
%     input is always 1 because makbalancedz assumes that samples or trials
%     are indexed across rows.
% 
% 
% Output
%   
%   zb - Balanced z-scores after transforming the values in x. zb will have
%     the same type as x.
%   
%   mbal - The balanced mean for each grouping of rows according to g
%   
%   sbal - The balanced standard deviation for each grouping of rows
%     according to g
% 
% 
% Balanced z-scoring is a step towards computing corrected choice or detect
% probabilities that would otherwise be under-estimated. The issue arrises
% when there is an inequal number of trials for the two types of
% behavioural response. This will tilt the estimated mean and standard
% deviation of a conventional z-scoring towards one sub-set of trials or
% another, blurring the threshold between the distribution of responses on
% each trial type. The trick is to pretend that there are actually an even
% number of trials and to combine the mean and standard deviation of the
% two sub-sets accordingly.
% 
% If behavioural responses of types 1 and 2 produce sub-sets of n1 and n2
% trials that each have estimated means of m1 and m2 and standard
% deviations of s1 and s2, then the composite mean mcom and standard
% deviation scom are found by computing:
%   
%   mcom = ( n1 * m1  +  n2 * m2 )  /  ( n1 + n2 )
%   
%   scom = sqrt (  A  +  B  )
%   
%     where  A = ( n1 * s1 ^ 2  +  n2 * s2 ^ 2 )  /  ( n1 + n2 )
%            
%            B = n1 * n2 * ( m1 - m2 ) ^ 2  /  ( n1 + n2 ) ^ 2
% 
% To compute the balanced mean mbal and standard deviation sbal, let 
% n1 = n2 in the equations above to obtain:
% 
%   mbal = ( m1  +  m2 )  /  2
%   
%   sbal = sqrt (  A  +  B  )
%   
%     where  A = ( s1 ^ 2  +  s2 ^ 2 )  /  2
%            
%            B = ( m1 - m2 ) ^ 2  /  4
% 
% If the balanced mean and standard deviation is used to compute balanced
% z-scores for trials grouped by stimulus condition then all balanced z-
% scores can be combined into a single set to compute CP or DP.
% 
% 
% Reference
% 
%   Kang, I. and J. H. Maunsell (2012). "Potential confounds in estimating
%     trial-to-trial correlations between neuronal response and behavior
%     using choice probabilities." J Neurophysiol 108(12): 3403-3415.
% 
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Maximum number of inputs / outputs
   narginchk (  3  ,  5  )
  nargoutchk (  0  ,  3  )
  
  % Size of x
  sx = size (  x  ) ;
  
  % We will allow makfun to do the work of input checking of x and g. While
  % mean and std can check outtype and nanflag. Here , we will check if b
  % is a logical vector with a length equal to the size of dim 1 in x.
  if  ~ isvector (  b  )  ||  ~ islogical (  b  )  ||  ...
      length (  b  )  ~=  sx( 1 )
    
    error (  'MAK:makbalancedz:b'  ,  [ 'makbalancedz: b must be a ' , ...
      'logical vector with length equal to size of x dim 1' ]  )
    
  % And check that g has as many rows as x
  elseif  size (  g  ,  1  )  ~=  sx( 1 )
    
    error (  'MAK:makbalancedz:g'  ,  [ 'makbalancedz: g must have ' , ...
      'as many rows as x' ]  )
    
  end % check input
  
  % Guarantee that b is a column vector
  if  ~ iscolumn (  b  )  ,  b = b( : ) ;  end
  
  
  %%% Balanced z-scoring %%%
  
  % Determine sub-sets of trials under all grouping conditions
  [ G , s ] = makfun (  [  g  ,  b  ]  ) ;
  
  % Compute mean and standard deviation of each sub-group
  m = makfun (  @( x ) mean( x , 1 , varargin{ : } )  ,  x  ,  G  ,  s  ) ;
  sd = ...
   makfun (  @( x ) std( x , 0 , 1 , varargin{ : } )  ,  x  ,  G  ,  s  ) ;
 
  % Let's guarantee that m and sd have a single row. This will be the case
  % if the size of x in dimensions 2 and higher are greater than 1. But it
  % will not be if x is a vector.
   m = reshape (   m  ,  [ 1 , sx( 2 : end ) , s ]  ) ;
  sd = reshape (  sd  ,  [ 1 , sx( 2 : end ) , s ]  ) ;
 
  % For flexible indexing , return all elements along all dimensions up to
  % but not including the balancing factor
  I = repmat (  { ':' }  ,  1  ,  ndims( m ) - 1  ) ;
  
  % Pull out m1, m2, s1, s2 for convenience
  m1 = m( I{ : } , 1 ) ;
  m2 = m( I{ : } , 2 ) ;
  
  s1 = sd( I{ : } , 1 ) ;
  s2 = sd( I{ : } , 2 ) ;
 
  % Compute the balanced mean ...
  mbal = (  m1  +  m2  )  ./  2 ;
  
  % ... and balanced standard deviations
  A = (  s1 .^ 2  +  s2 .^ 2  )  ./  2 ;
  
  B = ( m1 - m2 ) .^ 2  ./  4 ;
  
  sbal = sqrt (  A  +  B  ) ;
  
  % If all entries in a column of x are zeros then the balanced standard
  % deviation will be zero. And yet, the resulting balanced z-scores will
  % be NaN following division by zero. Therefore, any balanced standard
  % deviations equal to zero will be replaced by 1's so that zeros are
  % returned as balanced z-scores.
  sbal(  sbal  ==  0  ) = 1 ;
  
  % We now identify all trials in each sub-group of rows identified by g ,
  % but not b
  G = G( : , 1 : end / 2 )  |  G( : , end / 2 + 1 : end ) ;
  
  % Allocate output
  zb = zeros (  sx  ,  'like'  ,  x  ) ;
  
  % Reshape I to access all remaining columns in x
  I = repmat (  { ':' }  ,  1  ,  ndims( x ) - 1  ) ;
  
  % Sub-sets of trials
  for  i = 1 : prod (  s( 1 : end - 1 )  )
    
    % Point to sub-set of trials for this particular set of grouping factor
    % values
    j = G( : , i ) ;
    
    % Subtract mean
    zb( j , I{ : } ) = bsxfun (  @minus    ,  ...
       x( j , I{ : } )  ,  mbal( 1 , I{ : } , i )  ) ;
    
    % Divide by standard deviation
    zb( j , I{ : } ) = bsxfun (  @rdivide  ,  ...
      zb( j , I{ : } )  ,  sbal( 1 , I{ : } , i )  ) ;
    
  end % trial sub-sets
  
  % mbal requested , squeeze it down
  if  1  <  nargout  ,  mbal = squeeze (  mbal  ) ;  end
  
  % sbal requested , squeeze it down
  if  2  <  nargout  ,  sbal = squeeze (  sbal  ) ;  end
  
  
end % makbalancedz

