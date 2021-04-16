
function  [ m , e , varargout ] = makbindvar (  x  ,  y  ,  varargin  )
% 
% [ m , e ] = makbindvar (  x  ,  y  )
% [ m , e ] = makbindvar (  x  ,  y  ,  nbins  )
% [ m , e ] = makbindvar (  x  ,  y  ,  edges  )
% [ m , e , b , C ] = makbindvar (  x  ,  y  ,  ...  )
% [ m , e ] = makbindvar (  ...  ,  'ValFun'  ,  vfun  )
% [ m , e ] = makbindvar (  ...  ,  'ErrFun'  ,  efun  )
% [ ... ] = makbindvar (  ...  ,  'Domain'  ,  dom  )
% 
% MET Analysis Kit. Bins the values of dependent variable y based on the
% values of independent variable x and returns the mean ( m ) and standard
% error of the mean ( e ) for each bin.
% 
% 
% Required input:
% 
% x - N element numeric vector. The independent variable used to define M
%   bins.
% 
% y - Numeric array with N rows and any size in higher dimensions. Rows of
%   y must be in register with elements of x. That is, the ith element of
%   x, x( i ), corresponds to the ith row of y, y( i , : , : , ... , : ).
%   This is the dependent variable that is binned according to the bin
%   assignment of x.
% 
%
% Optional input:
% 
% nbins - Number of bins to use. M == nbins.
% 
% edges - Vector of bin edges. Each element defines the left edge of each
% bin, except for the last element which defines the right edge of the
% final bin. The number of bins becomes numel( edges ) - 1.
% 
% If nbins/edges is ommited then the M bins are automatically chosen with
% uniform width bins that span the range of x in the chosen domain.
% 
% 
% Name/Value pairs:
% 
% 'ValFun' , mfun - Instead of returning the mean of each bin in m, use the
%   function provided by function handle mfun. Default @mean.
% 
% 'ValErr' , efun - Instead of returning the standard error of the mean in
%   e, use the function provided by the function handle efun. Default
%   @( x ) std( x ) ./ sqrt( size( x , 1 ) ).
% 
% 'Domain' , dom - Indicate the domain to use when defining bin edges.
%   Here, domain means the numeric range of x ('num', default) or the
%   percentile range of x ('per'). The way that nbins and edges are
%   interpreted changes according to domain:
%   
%   'num' - If bin edges are automatically detected or nbins is given then
%     M uniform width bins spanning the range of values in x are defined.
%     When edges are given then bin edges are taken to be values of x and
%     are directly compared to values of x. Thus, bins evenly span the
%     range of x but can result in large differences in the number of data
%     points per bin.
%   
%   'per' - Bin edges are defined using percentiles of x. In this case, if
%     automatic binning is used or nbins is given then M bins of uniform
%     width spanning from percentile 0 to 100. edges is interpreted as
%     being a set of percentile values spanning 0 to 100. Bin edges in
%     percentiles are then mapped to values in the numeric range of x, in
%     order to actually bin x. Depending on the distribution of x, this can
%     result in bins with variable widths but with an approximately even
%     number of data points.
% 
% 
% Output:
% 
% m - M row numeric array. Mean of each bin (row). Or the value of each bin
%   computed with function mfun. The size of dim 2 onwards matches the size
%   of y, assuming the mean is returned or that mfun operates within
%   columns (e.g. mean, median, std, var, etc.).
% 
% e - M row numeric array. Standard error of each bin. Or the value of
%   each bin computed with function efun.
% 
% b - M element numeric vector. Contains the centre of each bin in the
%   numeric domain of x.
% 
% C - M element cell array vector. Each element contains the logical index
%   vector identifying elements of x belonging to each bin. Apply C to y in
%   order to get values of y associated with each bin.
% 
% Written by Jackson Smith - Novembre 2018 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % Name value pair for cellfun suite of functions to return cell arrays
  UF = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Check Input %%%
  
  % Basic input output number check
   narginchk (  2  ,  9  )
  nargoutchk (  0  ,  4  )
  
  % Input parser object
  ipo = inputParser ;
  
  % Required input arguments
  addRequired (  ipo  ,  'x'  ,  @isnumeric  )
  addRequired (  ipo  ,  'y'  ,  @( y ) validateattributes( y , ...
    {'numeric'} , { 'nrows' , numel( x ) } , 2 )  )
  
  % Optional input arguments
  addOptional (  ipo ,  'bins' ,  [] ,  @( e ) validateattributes( e , ...
    { 'numeric' } , { 'vector' , 'increasing' } , 3 )  )

  % Name/Value parameters
  addParameter (  ipo  ,  'ValFun'  ,  @mean  ,  ...
    @( p ) isa( p , 'function_handle' )  )
  addParameter (  ipo  ,  'ErrFun'  ,  ...
    @( x ) std( x ) ./ sqrt( size( x , 1 ) )  ,  ...
      @( p ) isa( p , 'function_handle' )  )
  addParameter (  ipo  ,  'Domain'  ,  'num'  ,  ...
    @( p ) ischar( ...
      validatestring( p , { 'num' , 'per' } , 'makbindvar' ) )  )
    
	% Parse input arguments
  parse (  ipo  ,  x  ,  y  ,  varargin{ : }  )
  
  % Get mean and error functions
  mfun = ipo.Results.ValFun ;
  efun = ipo.Results.ErrFun ;
  
  
  %%% Bin values of x %%%
  
  % Default nbins is empty
  nbins = [] ;  edges = [] ;
  
  % Number of bins or list of edges given?
  if  isempty (  ipo.Results.bins  )
    
    % No. histcounts optional input is empty
    hcin = { } ;
    
  % nbins given
  elseif  isscalar (  ipo.Results.bins  )
    
    % Assign bin number, empty edges
    nbins = ipo.Results.bins ;
    
    % Can't be zero or less, can't have fractional component
    if  nbins  <=  0  ||  mod (  nbins  ,  1  )
      error (  'MAK:makbindvar:nbins'  ,  ...
        'makbindvar: nbins must be a positive scalar integer'  )
    end
    
    % histcounts optional input is the number of bins
    hcin = {  nbins  } ;
    
  % Edges given
  else
    
    % Assign edges and number of bins
    edges = ipo.Results.bins ;
    
  end % nbins/edges
  
  % Use percentile domain
  if  strcmp (  ipo.Results.Domain  ,  'per'  )
    
    % No edges defined yet
    if  isempty (  edges  )
      
      % No number of bins given
      if  isempty (  nbins  )
        
        % Automatically find number of bins in numeric range of x
        [ ~ , edges ] = histcounts (  x  ) ;

        % Get number of bins
        nbins = numel (  edges  )  -  1 ;
      
      end % number of bins
      
      % Create bins in percentile domain
      edges = 0 : 100 / nbins : 100 ;
      edges( end ) = 100 ;
      
    % Make sure that edges are not outside of range 0 to 100
    elseif  edges( 1 )  <  0  ||  100  <  edges( end )
        
      error (  'MAK:makbindvar:edges_per'  ,  ...
      'makbindvar: percentile edges must be in range 0 to 100'  )
      
    end % edges
    
    % Map from percentiles to numeric values of x
    edges = prctile (  x  ,  edges  ) ;
    
  end % perc domain
    
  % Edges defined
  if  ~ isempty (  edges  )
    
    % histcounts optional input is the set of bin edges
    hcin = {  edges  } ;
    
  end % edges
  
  % Bin x , getting edge set and bin assignment of all values in x
  [ ~ , edges , xb ] = histcounts (  x  ,  hcin{ : }  ) ;
  
  % Assign number of bins if unknown
  if  isempty (  nbins  )  ,  nbins = numel (  edges  )  -  1 ;  end
  
  % Identify rows of y that go with each bin of x
  C = arrayfun (  @( i ) xb == i  ,  1 : nbins  ,  UF{ : }  )' ;
  
  
  %%% Evaluate binned y %%%
  
  % We need to expand across all other dimensions of y , so get a colon
  % operator for each one
  cop = repmat (  { ':' }  ,  1  ,  ndims( y ) - 1  ) ;
  
  % Mean/value function
  m = cellfun (  @( i ) mfun(  y( i , cop{ : } )  )  ,  C  ,  UF{ : }  ) ;
  
  % Error function
  e = cellfun (  @( i ) efun(  y( i , cop{ : } )  )  ,  C  ,  UF{ : }  ) ;
  
  % Collapse into arrays
  m = cell2mat (  m  ) ;
  e = cell2mat (  e  ) ;
  
  
  %%% Return output %%%
  
  % Bin centres requested
  if  2  <  nargout
    
    % Centres
    b = (  edges( 1 : end - 1 )  +  edges( 2 : end )  )  /  2 ;
    
    % Return
    varargout{ 1 } = b( : ) ;
    
  end % bin centres
  
  % Bin assignments requested
  if  3  <  nargout
    
    % Return
    varargout{ 2 } = C ;
    
  end % bin assign
  
  
end % makbindvar

