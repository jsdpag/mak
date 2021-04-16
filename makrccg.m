
function  [ rccg , varargout ] = makrccg (  varargin  )
% 
% [ rccg , lags ] = makrccg (  w  ,  A  ,  B  ) - Implements Bair et al.
%   2001's r_ccg metric. This is a cross-correlation based measure of the
%   spike train correlation between two sets of spike trains, A and B ; it
%   can be used to identify the time scale of correlation between the two
%   sets, to millisecond precision. Each set is a cell array vector with
%   the same number of elements, one per trial. Each element of A or B will
%   contain a row vector of spike times, in seconds. w is a 2-element
%   vector defining an analysis window such that w( 1 ) is the starting
%   time of the window and w( 2 ) is the ending time, both in seconds.
%   Spikes from A and B that fall within this window are analysed.
%   
%   rccg( t + 1 ) is the r_ccg between A and B integrated to a lag of t
%   milliseconds, for any value of t from 0 to ceil( diff( w ) / 1000 ) - 1
%   i.e. the number of milliseconds minus 1 in the analysis window w. lags
%   is a vector the size of rccg containing the number of milliseconds of
%   lag that was integrated to compute each rccg value. That is, lags( i )
%   equals t when i equals t + 1.
% 
% [ rccg , lags ] = makrccg (  w  ,  C  ) - When C is an N by M cell array
%   of spike times. Rows are indexed over N trials while columns are
%   indexed over M sets of spike trains ; each column may represent a
%   unique spike cluster. r_ccg is computed between all M ^ 2 pairs of
%   clusters. rccg becomes a L x M x M 3D matrix where lags are indexed
%   across rows, and pairs of spike clusters are indexed across the second
%   and third dimension. Thus, rccg( : , i , j ) is the full r_ccg curve
%   between clusters i and j. rccg( : , i , i ) contains the shift-
%   corrected auto-correlations at all integration lags ; these are not
%   normalised and hence have a unit of spikes-per-millisecond squared ; in
%   the limit, it is the auto-covariance and can be used to recover the
%   cross-covariance from any two units at any given integration lag. lags
%   is an L x 1 vector of milliseconds of lag, in register with the rows of
%   rccg.
% 
% [ rccg , lags , nscx ] = makrccg (  w  ,  ...  ) - For any previous form
%   of makrccg, an optional third output argument can be returned. This
%   will be the shift-corrected cross-correlations. That is, it returns the
%   curves that are integrated then normalised to obtain rccg. nscx will
%   have the same size as rccg, except in dimension 1 which will have
%   length ( 2L - 1 ) spanning cross-correlation lags -L + 1 : +L - 1. Lags
%   for nscx can be obtained from lags by:  [  - lags( end : -1 : 1 ) ;
%                                              + lags( 2 : end )      ]
% 
% [ rccg , lags , P , X ] = makrccg (  w  ,  ...  ) for either of the two
%   forms above, makrccg can return optional intermediate data that is
%   required to compute r_ccg. P is an N x Q x S matrix of peri-stimulus
%   time histograms using millisecond-wide time bins. It indexes rows over
%   over N trials, columns over the Q bins that fit within window w
%   ( rounding up ), and layers over S spike clusters. X returns all auto-
%   and cross-correlations between all spike clusters for each trial. It is
%   an N x L x S matrix indexing N trials over rows, L millisecond lags
%   over columns, and S spike clusters over layers. P and X are used in the
%   following form of makrccg.
%   
%   NOTE! For this form to be activated both P and X must be requested in
%   the left-hand side of output arguments. If only P is requested then
%   normalised, shift-corrected cross correlograms are returned instead.
%   
%   NOTE! All forms of makrccg (  w  ,  ...  ) function call require the
%   Parallel processing toolbox.
%   
%   NOTE! Returning P and X consumes an exponentially increasing amount of
%   memory with each additional spike cluster, due to the expansion of
%   cluster pairings. If P and X are not requested as output arguments in
%   any makrccg (  w  ,  ...  ) form of the function call then makrccg will
%   accumulate data in such a way as to minimise memory use, which will be
%   necessary for large data sets. This means that there will be some small
%   difference in the numerical result, due to finite-precision roundoff
%   error.
% 
% rccg = makrccg (  P  ,  X  ) computes the r_ccg for all spike cluster
%   pairs. This form is useful when r_ccg must be computed for different
%   subsets of trials, or when bootstrapping must be performed. For
%   instance, one may call bootci (  2e3  ,  @makrccg  ,  P  ,  X  ) to get
%   the 95% BCA bootstrap confidence intervals around the r_ccg for every
%   spike cluster pair.
% 
% 
% Reference:
% 
%   Bair W, Zohary E, Newsome WT. 2001. J Neurosci. 21(5):1676-97.
% 
% 
% Part of MET Analysis Kit (MAK)
% 
% Written by Jackson Smith - March 2018 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % makrccg form codes
  WINARG = 1 ;
  PXARGS = 2 ;
  
  % Exact number of output arguments required to trigger return of nscx
  RETNSCX = 3 ;
  
  % UniformOutput false input argument
  UF = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Check input %%%
  
  % Check limits on number of args
  narginchk  (  2  ,  3  )
  nargoutchk (  0  ,  4  )
  
  % Get the form of the function call. If the first input argument is a
  % two-element numeric vector , then this is one of the first two forms of
  % makrccg.
  if  isvector (  varargin{ 1 }  )  &&  ismatrix (  varargin{ 1 }  )  &&...
        isnumeric (  varargin{ 1 }  )  &&  numel (  varargin{ 1 }  )  ==  2
          
    % Set makrccg form flag
    frmflg = WINARG ;
    
    % Get w
    w = varargin{ 1 } ;
    
    % Window width in seconds
    Nlags = diff (  w  ) ;
    
    % Second element must be greater than the first by at least one
    % millisecond
    if  w( 2 )  <=  w( 1 )
      
      error (  'MAK:makrccg:worder'  ,  ...
        'makrccg: w must be a two element vector of increasing values'  )
      
    % Minimum window width
    elseif  Nlags  <  0.002
      
      error (  'MAK:makrccg:wwidth'  ,  ...
        'makrccg: w must define a window at least 2 milliseconds wide'  )
      
    end % w too narrow
    
    % Number of lags in milliseconds
    Nlags = ceil (  Nlags  /  0.001  )  -  1 ;
    
    % Get spike cluster trains from C
    if  nargin  ==  2
      
      % Assign name
      C = varargin{ 2 } ;
      
      % Must have no more than 2 dimensions with at least 2 columns
      if  ~ iscell (  C  )  ||  isempty (  C  )

        error (  'MAK:makrccg:Ccell'  ,  ...
          'makrccg: C must be a non-empty cell array'  )

      elseif  ~ ismatrix (  C  )  ||  size (  C  ,  2  )  <  2

        error (  'MAK:makrccg:Csize'  ,  [ 'makrccg: C must have at ' , ...
          'least 2 columns and at most 2 dimensions' ]  )

      end % check C
      
    % Get spike cluster trains from A and B
    elseif  nargin  ==  3
      
      % Assign names
      A = varargin{ 2 } ;
      B = varargin{ 3 } ;
      
      % Must both be cell arrays
      if  ~ iscell (  A  )  ||  ~ iscell (  B  )  ||  ...
          isempty (  A  )  ||  isempty (  B  )
        
        error (  'MAK:makrccg:ABcell'  ,  ...
          'makrccg: A and B must be non-empty cell arrays'  )
        
      % Must be vectors
      elseif  ~ isvector (  A  )  ||  ~ isvector (  B  )
        
        error (  'MAK:makrccg:ABvector'  ,  ...
          'makrccg: A and B must be vectors'  )
        
      % Must have same number of elements
      elseif  numel (  A  )  ~=  numel (  B  )
        
        error (  'MAK:makrccg:ABnumel'  ,  ...
          'makrccg: A and B must have the same number of elements'  )
        
      end % check A B
      
      % Pack into C
      C = [  A( : )  ,  B( : )  ] ;
      
    end % spk clust trains
    
    % Check that all spike trains are numeric row vectors or empty
    c = cellfun (  ...
      @( c ) isnumeric( c ) && ( isrow( c ) || isempty( c ) )  ,  C  ) ;
    
    if  ~ all (  c( : )  )
      
      error (  'MAK:makrccg:spkrowvect'  ,  ...
        'makrccg: all spike trains must numeric row vectors or empties'  )
      
    end
    
    % Size of C
    csize = size (  C  ) ;
    
    % Number of trials , number of spike clusters
    Ntrials = csize (  1  ) ;
    Nclusts = csize (  2  ) ;
    
    
  % First input argument is a 3D matrix , this is the P X form of makrccg
  elseif  isnumeric (  varargin{ 1 }  )  &&  ...
      ndims (  varargin{ 1 }  )  ==  3  &&  nargin  ==  2
    
    % Too many output arguments
    if  1  <  nargout
      
      error (  'MAK:makrccg:argoutnum'  ,  ...
        'makrccg: only rccg returned for P and X input argumets'  )
      
    end
    
    % Set makrccg form flag
    frmflg = PXARGS ;
    
    % Get P and X
    P = varargin{ 1 } ;  X = varargin{ 2 } ;
    
    % Size of P
    psize = size (  P  ,  2  ) ;
    
    % Number of lags and spike clusters
    Nlags   = psize( 2 )  -  1 ;
    Nclusts = psize( 3 ) ;
    
    % P has at least two time bins
    if  psize( 2 )  <  2
      
      error (  'MAK:makrccg:Pbins'  ,  ...
        'makrccg: P must have at least 2 time bins'  )
    
    % Check X is numeric
    elseif  ~ isnumeric (  X  )
      
      error (  'MAK:makrccg:Xnumeric'  ,  'makrccg: X must be numeric'  )
      
    % Check size of X
    elseif  ~ all(  size (  X  )  ==  ...
        [ psize( 1 ) , 2 * Nlags + 1 , psize( 3 ) ]  )
      
      error (  'MAK:makrccg:Xsize'  ,  'makrccg: X is the wrong size'  )
      
    end % check X
    
    
  % First argument is unrecognised
  else
    
    error (  'MAK:makrccg:input'  ,  'makrccg: invalid input arguments'  )
    
  end % Get form of function call
  
  
  %%% Pre-processing %%%
  
  % w input arg form of makrccg , thus we need to do pre-processing
  if  frmflg  ==  WINARG
    
    % lags output argument requested
    if  1  <  nargout  ,  varargout{ 1 } = ( 0 : Nlags )' ;  end
    
    % We need to build P by binning spike times using millisecond-width
    % time bins
    edges = ( 0 : Nlags + 1 ) / 1e3  +  w( 1 ) ;
    
    % At P and X are requested as output , compute P and X for each trial
    % and use up lots of memory
    if  RETNSCX  <  nargout

      % Bin all spikes , P is Nclusts x Ntrials and contains column vects
      P = cellfun (  @( c ) histcounts( c , edges )'  ,  C  , UF{ : }  )' ;

      % Reshape cell array so that clusters span columns and trials span
      % layers
      P = reshape (  P  ,  1  ,  Nclusts  ,  Ntrials  ) ;

      % And collapse into a single matrix of time bins x clusters x trials
      P = cell2mat (  P  ) ;

      % Now we need to compute auto and cross correlations for each trial
      parfor  i = 1 : Ntrials

        X( : , : , i ) = xcorr (  P( : , : , i )  ) ;

      end % trials

      % Now permute P and X so that trials span rows
      P = permute (  P  ,  [ 3 , 1 , 2 ]  ) ;
      X = permute (  X  ,  [ 3 , 1 , 2 ]  ) ;

      % Return P and X
      varargout{ 2 } = P ;

      if  3  <  nargout
        varargout{ 3 } = X ;
      end
    
    % Neither P nor X are requested as output , so compute them in a way
    % that minimises memory usage
    else
      
      % Allocate P and X. P indexes bins over rows and spike clusters over
      % columns. X indexes lags over rows and spike cluster pairs over
      % columns i.e. cross- and auto-correlations.
      P = zeros (  numel( edges ) - 1  ,  Nclusts  ) ;
      X = zeros (  2 * ( Nlags + 1 ) - 1  ,  Nclusts ^ 2  ) ;
      
      % Trials
      parfor  i = 1 : Ntrials
        
        % Bin all spikes on this trial , p is 1 x Nclusts cell array column
        % vectors
        p = cellfun (  @( c ) histcounts( c , edges )'  ,  C( i , : )  ,...
          'UniformOutput'  ,  false  ) ;
        
        % Collapse p into a number of bins x spike clusters matrix of
        % spikes binned from this trial
        p = cell2mat (  p  ) ;
        
        % Accumulate spike count in each bin
        P = P  +  p ;
        
        % Accumulate cross-/auto-correlations
        X = X  +  xcorr (  p  ) ;
        
      end % trials
      
      % Average PSTH and cross-/auto-correlations
      P = P  ./  Ntrials ;
      X = X  ./  Ntrials ;
      
    end % build P and X
    
  end % w input arg form
  
  % If P and X contain data for each trial , either because they were
  % requested as output or because this is the makrccg (  P  ,  X  ) form
  % of the function call
  if  frmflg  ==  PXARGS  ||  RETNSCX  <  nargout
    
    % Compute the average PSTH and auto-/cross-correlation , over trials
    P = mean (  P  ,  1  ) ;
    X = mean (  X  ,  1  ) ;

    % Reshape so that bins/lags span rows and clusters span columns
    P = squeeze (  P  ) ;
    X = squeeze (  X  ) ;
    
  end % trial-wise P and X
  
  
  %%% Calculate r_ccg %%%
  
  % Subtract shift-predictors from all auto- and cross-correlations , at
  % every lag , we have not integrated yet
  X = X  -  xcorr (  P  ) ;
  
  % Get row indices of shift-corrected auto- and cross-correlations that
  % put each positive and negative lag pair in register with each other ,
  % ordered by absolute magnitude
  i = Nlags : -1 : 1 ;
  j = Nlags + 2 : 2 * Nlags + 1 ;
  
  % Sum values at each matching positive and negative lag i.e. find
  % C_a,b( -t )  +  C_a,b( +t ) , where C_a,b is the cross-correlation of
  % spike clusters a and b.
  A = [  X( Nlags + 1 , : )  ;  X( i , : ) + X( j , : )  ] ;
  
  % Cumulative sum across rows now provides integrated auto- and cross-
  % correlations for each possible lag
  A = cumsum (  A  ,  1  ) ;
  
  % Get the linear index of the lower-triangular portion of a Nclust x
  % Nclust matrix , excluding diagonal. This is a bit misleading because we
  % will use it to access the upper-triangular portion of the "NxN" matrix
  % of cross correlations returned by xcorr. Why? See doc xcorr:
  % 
  %   if S is a three-channel signal, S = ( x1 , x2 , x3 ) then the result
  %   of R = xcorr (  S  ) is organised as R = ( R11 , R12 , R13 , R21 ,
  %   R22 , R23 , R31 , R32 , R33 ).
  % 
  % But Matlab uses a column-major organisation , meaning that find looks
  % along columns first. If Nclusts is 3 then the linear index of the
  % output of tril is [ 2 , 3 , 6 ] , which accesses R12, R13, and R23 from
  % the above example i.e. the upper-triangular portion of a 3 x 3 matrix
  % with arrangement [ R11, R12, R13 ; R21, R22, R23 ; R31, R32, R33 ].
  i = find ( tril(  true( Nclusts )  ,  - 1  ) ) ;
  
  % Get the sub-script row and column indices for each linear index. Now we
  % have vectors i , a , b that are all the same size. i finds the cross-
  % correlation integral for each unique pair of spike clusters, while a
  % and b identify each spike cluster in the pair. See above for the reason
  % that we assign row subscripts to b and column subscripts to a.
  [ b , a ] = ind2sub (  [ Nclusts , Nclusts ]  ,  i  ) ;
  
  % However , we must then convert a and b into equivalent linear column
  % indices for A. That is, the linear indices of the diagonals. This
  % allows us to get the auto-correlations. We keep the subscript indices
  % in order to perform the transpose on rccg.
  ad = sub2ind (  [ Nclusts , Nclusts ]  ,  a  ,  a  ) ;
  bd = sub2ind (  [ Nclusts , Nclusts ]  ,  b  ,  b  ) ;
  
  % Allocate rccg memory
  rccg = zeros ( size(  A  ) ) ;
  
  % Compute r_ccg for each unique pair of spike clusters. Just to make it
  % interesting, while i accesses the upper-triangular part of the
  % correlation matrix from xcorr, it then places the result in the lower-
  % triangular part of rccg.
  rccg( : , i ) = A( : , i )  ./  sqrt (  A( : , ad )  .*  A( : , bd )  ) ;
  
  % Now find the linear index of the upper-triangular portion of a Nclust x
  % Nclust matrix
  j = sub2ind (  [ Nclusts , Nclusts ]  ,  a  ,  b  ) ;
  
  % And copy all r_ccg values to the other side of the correlation matrix ,
  % sort of like a transposition-copy
  rccg( : , j ) = rccg( : , i ) ;
  
  % Locate all rccg auto-correlations
  i = logical(  eye( Nclusts )  ) ;
  
  % Return shift-corrected auto-correlations at all integration lags
  rccg( : , i ) = A( : , i ) ;
  
  % nscx requested , there's some extra work to do
  if  nargout  ==  RETNSCX
    
    % First , compute the average firing rate of all cells within the given
    % time window
    parfor  i = 1 : Nclusts
      
      G( 1 , i ) = mean (  ...
        cellfun (  @( c )  sum (  w( 1 ) <= c  &  c <= w( 2 )  )  ,  ...
          C( : , i )  )  )  /  diff (  w  ) ; %#ok
      
    end % average rate
    
    % Now compute the geometric mean firing rate of each pair ,
    % transposition puts linear indexing of geo rate in register with
    % columns of X
    G = sqrt (  bsxfun(  @times  ,  G'  ,  G  )  )' ;
    
    % Triangular shape function removes the CCG artefact caused by zero
    % padding
    T = ( Nlags + 1 )  -  [ Nlags : -1 : 0 , 1 : Nlags ]' ;
    
    % Normalise the shift-corrected cross-correlations by the geo mean AND
    % triangular function
    X( : ) = X  ./  bsxfun (  @times  ,  T  ,  G( : )'  ) ;
    
  end % nscx request
  
  % Special case: ( w , A , B ) input argument list
  if  frmflg == WINARG  &&  2 < nargin
    
    % return cross-correlation rccg vector
    rccg = rccg (  :  ,  2  ) ;
    
    % nscx requested , return vector
    if  nargout  ==  RETNSCX
      varargout{ RETNSCX - 1 } = X (  :  ,  2  ) ;
    end
    
  % Standard case
  else
    
    % Reshape into lags x clusters x clusters matrix
    rccg = reshape (  rccg  ,  Nlags + 1  ,  Nclusts  ,  Nclusts  ) ;
    
    % nscx requested , reshape similarly to rccg. The 2nd and 3rd
    % dimensions are transposed so that X( : , m , n ) provides the cross-
    % correlation of spike cluster m and spike cluster n where n is shifted
    % by each lag ; this affects how the cross-correlogram is relfected
    % across lag zero as X( : , m , n ) == X( end : -1 : 1 , n , m ).
    if  nargout  ==  RETNSCX
      X = reshape (  X  ,  2 * Nlags + 1  ,  Nclusts  ,  Nclusts  ) ;
      X = permute (  X  ,  [ 1 , 3 , 2 ]  ) ;
      varargout{ RETNSCX - 1 } = X ;
    end
    
  end % cases
  
end % makrccg

