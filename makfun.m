
function  varargout = makfun (  varargin  )
% 
% [ y , n , u ] = makfun (  h  ,  x  ,  v  )  groups rows in x by unique
%   values in the grouping variable v, and then applies function handle h
%   to each grouping of x. x can be a numeric, char, logical, cell, or
%   struct array ; and v must be numeric, char, or logical. x may be N-
%   dimensional. v must be a vector whose length equals the number of rows
%   in x, thus assigning a group to each row. The function pointed to by h
%   must return a numeric array that can be M-dimensional. The results are
%   collected in y, which has M + 1 dimensions, where the final dimension
%   is indexed over the ascending unique values of v. In other words, it is
%   indexed over groupings of x. n is a vector with the number of values
%   per grouping of x. If x is empty then y is empty and n has only zeros.
%   u is the vector of unique values in v, sorted ascending.
% 
% [ y , n , U ] = makfun (  h  ,  x  ,  v_1  ,  v_2  , ... v_K  )  groups
%   values in x by each combination of unique values in grouping variables
%   v_1 through v_K, and then applies the function in handle h to each
%   group. All grouping variables must be vectors with length equal to the
%   number of rows in x. y becomes an M + K dimensional array there the M +
%   ith dimension is indexed over ascending unique values of grouping
%   variable v_i. n becomes a K dimensional array with size [ u_1 , u_2 ,
%   ... u_k ], where u_i is the number of unique values in v_i. U is a cell
%   array with the set of unique values for each grouping variable.
% 
% [ y , n , U ] = makfun (  h  ,  x  ,  V  )  is an alternative form of
%   makfun( h , x , v_1 , v_2 , ... v_K ) where V is a matrix such that
%   V = [ v_1 , v_2 , ... v_K ], and all v_i are column vectors.
% 
% [ G , s , n , U ] = makfun (  v_1  , v_2  , ... v_K  )
% [ G , s , n , U ] = makfun (  V  )  When h and x are excluded and only
%   the grouping variables are provided, then makfun returns grouping
%   information. n and U are the same as above. In reverse order, s is a
%   length K vector with the number of unique values in each grouping
%   variable ; hence s( i ) is the number of unique values in v_i or 
%   V( : , i ). It follows that numel( U{ i } ) == s( i ). All grouping
%   variables have the same number of elements, let this number be P. G
%   then becomes a logical index matrix of size P x ( s( 1 ) * s( 2 ) * ...
%   s( K ) ) i.e. P x prod( s ), where the ith column of G is the subset of
%   elements that are grouped by the ith combination of unique grouping
%   variable values. One obtains grouping information when the same
%   grouping will be applied in repeated calls to makfun with different
%   values for h or x. This saves time by allowing the grouping information
%   to be calculated just once, rather than for every call to makfun.
% 
% y = makfun (  h  ,  x  ,  G  ,  s  )  does the same thing as a
%   call to makfun( h , x , v_1 , v_2 , ... v_K ) or makfun( h , x , V ). x
%   and G must have the same number of rows.
% 
% [ Y , ... ] = makfun (  H  ,  x  ,  ...  )  each function handle in cell
%   array H is applied to every grouping of x. Y is a cell array with the
%   same size as H, where Y{ i } contains the output of function H{ i }. n
%   remains a length K vector with the number of values per grouping. If x
%   is empty in this case then Y is a cell array of empty numeric matrices.
% 
% ... = makfun (  h  ,  X  ,  ...  )  if function handle h requires N input
%   arguments then these can be provided in cell array X. X must be an
%   N-element row vector ; if it has more than one row then it is
%   considered to be the function data , and its rows will be grouped into
%   subsets and passed directly to h. Each 
%   element of X must have the same number of rows, but can be of any valid
%   type for x ( see above ). The same grouping is applied to all elements
%   of X before passing them as arguments to h. Any grouping variables
%   v_<i>, grouping  variable matrix V, or grouping information matrix G
%   must have the same number of rows as the matrices in X. If h is a cell
%   array of function handles, then each function must accept at least N
%   input arguments.
% 
% Part of MET Analysis Kit (MAK)
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % Drops input argument pair into array- and cellfun , telling them to
  % return a cell array. Pass in argument as comma-separated list by
  % putting UF{ : } as the last input argument to either function.
  UF = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Check input %%%
  
  % Input argument variables , initialised to empty
  h = [] ;  x = [] ;  rx = [] ;  G = [] ;  s = [] ;
  
  % No input arguments
  if  nargin  ==  0
    error (  'MAK:makfun:noargsin'  ,  'makfun: input argument required'  )
  end
  
  
  %-- Function output form of makfun --%
  
  % Determine form of the function call. Do we return function output or
  % grouping information? We know it is the former if argument 1 is a
  % function handle or a cell array of function handles.
  if  isfun (  varargin{ 1 }  )
    
    % There must be at least 3 inputs and at most three outputs
    narginchk  (  3  ,  Inf  )
    nargoutchk (  0  ,    3  )
    
    % Get function handle(s) , and putative function data
    h = varargin{ 1 } ;
    x = varargin{ 2 } ;
    
    % If h is a function handle then pack it into a cell array so that we
    % can write generic code that covers all situations
    if  ~ iscell (  h  )  ,  h = { h } ;  end
    
    % For the same reason , check whether x is a cell array row vector of
    % multiple input arguments to h. If it is not, i.e. it is function data
    % for a single input argument, then pack in a container cell array.
    % Again, for generic code.
    if  ~ (  iscell (  x  )  &&  isrow (  x  )  )  ,  x = { x } ;  end
    
    % All elements of x must be numeric, char, logical, cell, or struct
    if  ~ all ( cellfun(  @chxtype  ,  x  ) )
      error (  'MAK:makfun:xtype'  ,  ...
        'makfun: x must be numeric, char, logical, cell, or struct'  )
    end
    
    % Get number of rows from the first set of function data in x
    rx = size (  x{ 1 }  ,  1  ) ;
    
    % Make sure that all function data has the same number of rows
    if  ~ all ( cellfun(  @( x ) size( x , 1 ) == rx  ,  x( 2 : end )  ) )
      error (  'MAK:makfun:xrows'  ,  ...
        'makfun: all elements of X must have the same number of rows'  )
    end
    
    % Is x followed by grouping information or grouping variables? We know
    % it is the former if the next input argument is a logical matrix.
    if  islogical (  varargin{ 3 }  )
      
      % In this case we require exactly four input arguments and at most
      % one output argument
      narginchk  (  4  ,  4  )
      nargoutchk (  0  ,  1  )
      
      % Fetch grouping information
      G = varargin{ 3 } ;
      s = varargin{ 4 } ;
      
      % Is s numeric?
      if  ~ isnumeric (  s  )
        
        error (  'MAK:makfun:snumeric'  ,  'makfun: s must be numeric'  )
        
      % Does s contain integers?
      elseif  any ( mod(  s  ,  1  ) )
        
        error (  'MAK:makfun:sints' ,  ...
          'makfun: s must only contain integers'  )
      
      % Convert s to double if it is not so already
      elseif  ~ isa (  s  ,  'double'  )
        
        s = double (  s  ) ;
        
      end % check s
      
      % Make s a row vector if it is not already
      if  ~ isrow (  s  )  ,  s = s( : )' ;  end
      
      % Size of G
      sg = size (  G  ) ;
      
      % Total number of groupings for all combinations of unique grouping
      % values across all grouping variables
      Ntotal = prod (  s  ) ;
      
      % Does G have as many rows as x does?
      if  sg( 1 )  ~=  rx
        
        error (  'MAK:makfun:G_rows' ,  ...
          'makfun: G must have as many rows as x'  )
        
      % Does G have as many columns as there are groupings?
      elseif  sg( 2 )  ~=  Ntotal
        
        error (  'MAK:makfun:G_rows' ,  ...
          'makfun: G must have as many columns as there are groupings'  )
        
      end % check G
      
    % Grouping variables expected
    else
      
      % Grouping variables start at third input argument
      i = 3 ;
      
    end % grouping info provided
    
  % Grouping information output
  else
    
    % Check number of output arguments
    nargoutchk (  0  ,  4  )
    
    % Grouping variables start at first argument
    i = 1 ;
    
  end % Check for [ y , n ] = makfun ( h , x , ... ) form
  
  
  %-- Grouping information --%
  
  % Grouping information was not provided , we expect grouping variables
  % from input argument i
  if  isempty (  G  )
    
    % Expand an index vector to get all remaining input arguments
    i = i : nargin ;
    
    % All remaining inputs must be numeric or char
    if  ~ all ( cellfun(  ...
          @( v ) isnumeric( v ) || ischar( v ) || islogical( v )  ,  ...
            varargin( i )  ) )
      
      error (  'MAK:makfun:gvarnumeric' ,  ...
          'makfun: grouping variables must be numeric or char'  )
        
    % Grouping variables are a set of vectors
    elseif  all ( cellfun(  @isvector  ,  varargin( i )  ) )
      
      % Check they all have the same length
      if  ~ all(  numel( varargin{ i( 1 ) } )  ==  ...
                cellfun( @numel , varargin( i( 2 : end ) ) )  )
        
        error (  'MAK:makfun:gvarlengths' ,  [ 'makfun: grouping ' , ...
        'variables vectors must all be the same length' ]  )
        
      end % same length
      
      % Get grouping variables in a cell array
      v = varargin( i ) ;
      
      % And guarantee that all are column vectors
      v = cellfun (  @( v ) v( : )  ,  v  ,  UF{ : }  ) ;
      
    % Grouping variables are columns of one matrix
    elseif  numel (  i  ) == 1  &&  ismatrix (  varargin{ i }  )  &&  ...
        ~ isvector (  varargin{ i }  )
      
      % Expand into a cell array of column vectors
      v = num2cell (  varargin{ i }  ,  1  ) ;
      
    % Invalid input
    else
      
      error (  'MAK:makfun:gvarinvalid' ,  [ 'makfun: grouping ' , ...
        'variables must be a set of vectors, each its own input ' , ...
        'argument, or a single matrix' ]  )
      
    end % check grouping variables
    
    % If this is the function output form of makfun then check that
    % grouping variables match the number of rows in x
    if  ~ isempty (  h  )  &&  numel (  v{ 1 }  )  ~=  rx
      
      error (  'MAK:makfun:gvarlengths' ,  [ 'makfun: grouping ' , ...
        'variables must have as many elements as there are rows in x' ]  )
      
    end % rows in x
    
    % Start computing grouping information. Start by getting unique set of
    % values from each grouping variable, sorted ascending.
    U = cellfun (  @unique  ,  v  ,  UF{ : }  ) ;
    
    % Count the number of unique values per grouping variable
    s = cellfun (  @numel  ,  U  ) ;
    
    % Total number of groupings for all combinations of unique grouping
    % values across all grouping variables
    Ntotal = prod (  s  ) ;
    
    % Locate each unique value in each grouping variable vector. Returns a
    % logical matrix for each grouping vector with columns indexed over the
    % vector's unique values. The ith column is a logical index that
    % locates all instances of the ith unique value in the grouping vector.
    % This is the first step to finding all possible groupings.
    G = cellfun (  @( v , u ) catmat( ...
      arrayfun(  @( u ) v == u  ,  u  ,  UF{ : }  ) )  ,  v  ,  U  ,  ...
        UF{ : }  ) ;
      
    % If there is only one grouping variable then we have already defined
    % all groups. Fetch the only matrix in G.
    if  isscalar (  v  )
      
      G = G{ 1 } ;
      
    % There are multiple grouping variables
    else

      % We might consider s to be the size of some N-dimensional array. Or,
      % we might consider there to be an N-dimensional space that spans the
      % unique values of all N grouping variables. Allocate a subscript
      % index cell vector.
      S = cell ( size(  G  ) ) ;

      % Get sub-script indices for the previously mentioned N-dimensional
      % array. This gives the unique value from each grouping variable for
      % each individual grouping. It may also be considered a set of N-
      % dimensional coordinates in the unique-value space.
      [  S{ : }  ] = ind2sub (  s  ,  1 : Ntotal  ) ;

      % Find logical index for each grouping. The ith column locates
      % elements that contain a specific combination of unique values
      % across all grouping variable vectors. This is the second and final
      % step of computing grouping information, building upon the first
      % assignment to G, above. We pass in the set of points in the N-
      % dimensional unique- value space. Each coordinate is received as a
      % set of scalar input arguments to the anonymous function of
      % arrayfun. Each scalar is passed in alongside the logical matrices
      % that locate unique values in each grouping variable ; each scalar
      % is a column index , as matrix columns are indexed over unique
      % values. The cellfun anonymous function retreives the column
      % locating one unique value for a given grouping variable, and
      % cellfun's output is a cell array of logical indices locating a set
      % of values across grouping variables. These are collapsed into a
      % matrix , and the logical index is returned that locates each
      % element in which every grouping variable has a specified value.
      G = arrayfun(  @( varargin ) all(  catmat(...
        cellfun(  @( g , c ) g( : , c ) ,  G ,  varargin ,  UF{ : }  )),...
          2 )  ,  S{ : }  ,  UF{ : }  ) ;

      % Collapse into a matrix
      G = catmat (  G  ) ;
    
    end % get grouping matrix
    
    % This is the grouping information form of makfun
    if  isempty (  h  )
      
      % Return G
      varargout{ 1 } = G ;
      
      % s requested
      if  1  <  nargout  ,  varargout{ 2 } = s ;  end
      
      % n requested
      if  2  <  nargout  ,  varargout{ 3 } = gnumel (  G  ,  s  ) ;  end
      
      % U requested
      if  3  <  nargout  ,  varargout{ 4 } = U ;  end
      
      % Done
      return
      
    end % grouping info form
    
  end % grouping vars
  
  
  %%% x is empty %%%
  
  if  isempty (  x  )
    
    % Return value depends on h. If h has a single function handle ...
    if  isscalar (  h  )
      
      % ... then return a plain empty
      varargout{ 1 } = [ ] ;
      
    % But if h has multiple handles ...
    else
      
      % ... then return a cell array of empties
      varargout{ 1 } = cell ( size(  h  ) ) ;
      
    end % check h
    
    % Return zeros vector , if requested
    if  1  <  nargout  ,  varargout{ 2 } = zeros ( [ s , 1 ] ) ;  end
    
    % Return unique grouping values , if requested
    if  2  <  nargout  ,  varargout{ 3 } = U ;  end
    
    % Done
    return
    
  end % empty x
  
  
  %%% Apply functions to grouped data %%%
  
  % Prepare a cell array of subscripts that will access all elements in all
  % dimensions greater than 1. This is done by passing in the colon
  % operator for the subscript of each dimension via comma-separated list.
  % Prepare a set of colon operators for each set of function data.
  S = cellfun (  @( x )  repmat(  { ':' } ,  ndims( x ) - 1 ,  1  )  ,  ...
    x  ,  UF{ : }  ) ;
  
  % Allocate output
  y = cell (  Ntotal  ,  numel( h )  ) ;
  
  % Groupings
  for  i = 1 : Ntotal
    
    % Group function data
    gx = cellfun(  @( x , S )  x(  G( : , i )  ,  S{ : }  )  ,  ...
      x  ,  S  ,  UF{ : }  ) ;
    
    % Apply functions to grouped data
    y( i , : ) = cellfun (  @( h ) h( gx{ : } )  ,  h  ,  UF{ : }  ) ;
    
  end % groupings
  
  % Get size vector for function output using placeholder values of 1 for
  % each non-singleton dimension. If the number of dimensions is 2 but the
  % dimension 2 has length 1 then the number of dimensions is revised to 1.
  % Assume that function output from the first grouping gives the size of
  % all other output.
  sy = cellfun (  @foutsize  ,  y( 1 , : )  ,  UF{ : }  ) ;
  
  % Collapse the output into numeric matrices , for each function. 
  y = cellfun (  ...
    @( y , sy ) makcell2mat(  reshape( y , [ sy , s , 1 ] )  )  ,  ...
      num2cell( y , 1 )  ,  sy  ,  UF{ : }  ) ;
  
  % Return function output for a single function handle ...
  if  isscalar (  h  )
    
    varargout{ 1 } = y{ 1 } ;
    
  % ... or for multiple function handles
  else
    
    varargout{ 1 } = y ;
    
  end % return y
  
  % Thanks to input checking we can only have more than one output argument
  % if grouping vectors were given as input arguments rather that pre-
  % computed grouping information. Return n if requested.
  if  1  <  nargout  ,  varargout{ 2 } = gnumel (  G  ,  s  ) ;  end
  
  % Return u if requested
  if  2  <  nargout  ,  varargout{ 3 } = U ;  end
  
  
end % makfun


%%% Subroutines %%%

% Checks if h is a function handle or a cell array of function handles
function  r = isfun (  h  )
  
  % If h is not a cell array then pack it into one
  if  ~ iscell (  h  )  ,  h = { h } ;  end
  
  % Check whether each element of cell array is a function handle
  c = cellfun (  @( h ) isa( h , 'function_handle' )  ,  h  ) ;
  
  % Return true if all elements held function handles
  r = all (  c  ) ;
  
end % isfun


% Checks that function data in x is of the correct type
function  r = chxtype (  x  )
  
  r = isnumeric( x )  ||  ischar( x )  ||  islogical( x )  ||  ...
    iscell( x )  ||  isstruct( x ) ;
  
end % chxtype


% Concatenate all elements of input cell array into an output matrix
function  m = catmat (  C  )
  
  m = [  C{ : }  ] ;
  
end % catmat


% Number of elements per grouping of x
function  n = gnumel (  G  ,  s  )
  
  % Count values per grouping
  n = sum (  G  ,  1  ) ;

  % Reshape according to number of grouping variables and number of
  % unique values within each
  n = reshape (  n  ,  [ s , 1 ]  ) ;
  
end % gnumel


% Function output size. An N-element row vector of ones is returned for the
% first N non-singleton dimensions of y.
function r = foutsize (  y  )
  
  % Get size of y
  s = size (  y  ) ;
  
  % Now find the last non-singleton dimension
  n = find (  1  <  s  ,  1  ,  'last'  ) ;
  
  % Return vector of placeholder 1's
  r = ones (  1  ,  n  ) ;
  
end % foutsize

