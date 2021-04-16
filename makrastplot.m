
function  [ varargout ] = makrastplot (  T  ,  varargin  )
%
% makrastplot (  T  )
% makrastplot (  ...  ,  Name  ,  Value  )
% h = makraster (  ...  )
% 
% MET Analysis Kit. Makes a raster plot in which sets of time stamps span
% the y-axis and time spans the x-axis. Each element of cell array T
% provides a vector of time stamps that occupy a single row. All time
% stamps are drawn using a single line object, the handle of which is
% optionally returned. By default, each time stamp is represented by a
% black dot '.' marker. Optional Name/Value pairs can be given to specify
% additional properties of the line object. The raster is added to the
% current axes, or a new figure and axes is created if none exist.
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% Default line properties %%%
  
  NameValue = {  'LineStyle'  ,  'none'  ,  'Marker'  ,  '.'  ,  ...
      'MarkerEdgeColor'  ,  'k'  } ;
    
  
  %%% Check input %%%
  
  % Number of arguments check
  narginchk (  1  ,  Inf  )
  nargoutchk (  0  ,  1  )
  
  % Check that T is a cell array
  if  ~ iscell (  T  )
    
    error (  'MAK:makrastplot:Tnotcell'  ,  ...
      'makrastplot: T must be a cell array'  )
    
  % All elements must be real-valued, numeric vectors with no Inf or NaN,
  % or empties
  elseif  ~ all ( cellfun(  @val  ,  T  ) )
    
    error (  'MAK:makrastplot:Tvalvect'  ,  [ 'makrastplot: ' , ...
      'All elements of T must be numeric, real-valued vectors with ' , ...
        'no Inf or NaN, or empty' ]  )
      
	% There are name/value pairs
  elseif  1  <  nargin
    
    % There is an even number of inputs , hence a complete number of pairs
    if  mod (  nargin - 1  ,  2  )
      
      error (  'MAK:makrastplot:NameValueUneven'  ,  ...
        'makrastplot: uneven number of property Names and Values'  )
    
    % First value of each pair is a string
    elseif  ~ all (  cellfun (  @isstring  ,  varargin( 1 : 2 : end )  )  )
      
      error (  'MAK:makrastplot:NameValueStrings'  ,  ...
        'makrastplot: Missing Name from Name/Value pairs'  )
      
    end
    
  end % check input
  
  
  %%% Preparation %%%
  
  % Number of time stamps per vector
  n = cellfun (  @numel  ,  T  ) ;
  
  % Guarantee that this is a row vector
  if  ~ isrow ( n )  ,  n = reshape (  n  ,  1  ,  numel( n )  ) ;  end
  
  % Get row number of each time stamp
  row = arrayfun (  @( r , n ) r * ones( 1 , n )  ,  ...
    1 : numel( T )  ,  n  ,  'UniformOutput'  ,  false  ) ;
  
  
  %%% Draw raster %%%
  
  % Current hold state
  h = ishold ;
  
  % Current plot hold is off , so refresh
  if  ~ h  ,  newplot  ,  end
  
  % Draw time stamps
  lh = line (  [ T{ : } ]  ,  [ row{ : } ]  ,  ...
    NameValue{ : }  ,  varargin{ : }  ) ;
  
  % Hold state is off , so tighten axes around tick marks
  if  ~ h  ,  axis tight  ,  end
  
  % Return line handle on request
  if  nargout  ,  varargout{ 1 } = lh ;  end
  
  
end % makrastplot


%%% Subroutines  %%%

% Checks if t is valid i.e. a numeric, real-numbered vector that contains
% no Inf or Nan values
function  x = val (  t  )
  
  % Valid vector
  vect = isnumeric (  t  )  &&  isreal (  t  )  &&  isvector (  t  )  &&...
    all ( isfinite(  t  ) ) ;
  
  % t is a valid vector or it is empty
  x = vect  ||  isempty (  t  ) ;
  
end % ists


% Checks that s is a char row vector
function  x = isstring (  s  )
  
  x = ischar (  s  )  &&  isrow (  s  ) ;
  
end % istring

