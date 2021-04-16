
function  [ m , s ] = makmedstd (  x  ,  varargin  )
% 
% [ m , s ] = makmedstd (  x  )
% ... = makmedstd (  x  ,  dim  )
% ... = makmedstd (  ...  ,  nanflag  )
% 
% MET Analysis Kit. Returns median of x in m and the standard deviation
% around the median in s. That is, the median is used to compute standard
% deviation from samples in x rather than the mean. Optional arguments dim
% and nanflag are the same as for the median function. dim says to compute
% the median and median-standard deviation across that dimension of x.
% nanflag can be 'omitnan' to ignore NaN values.
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  % Compute median of x
  m = median (  x  ,  varargin{ : }  ) ;
  
  % Get squared differences
  sd = bsxfun (  @minus  ,  x  ,  m  )  .^  2 ;
  
  % Number of values when dim given
  if  1  <  nargin  &&  isnumeric (  varargin{ 1 }  )
    
    % Grab dim
    dim = varargin( 1 ) ;
    
    % Number of values along that dimension
    n = size (  x  ,  dim{ 1 }  ) ;
    
  % When dim not given
  else
    
    dim = {} ;
    n = numel (  x  ) ;
    
  end % num vals
  
  % Remove NaN from sum if omitnan flag given
  if  1  <  nargin  &&  strcmp (  'omitnan'  ,  varargin{ end }  )
    
    % Count NaN values
    nn = sum (  isnan(  x  )  ,  dim{ : }  ) ;
    
    % Subtract them from number of values
    n = n  -  nn ;
    
  end % remove NaN

  % Square root of the normalised sum of squared differences
  s = sqrt (  sum (  sd  ,  varargin{ : }  )  /  ( n - 1 )  ) ;
  
  
end % makmedstd

