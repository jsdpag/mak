
function  varargout = maksttc (  varargin  )
% 
% [ sttc , dt ] = maksttc (  w  ,  maxdt  ,  A  ,  B  )
% [ sttc , dt ] = maksttc (  w  ,  maxdt  ,  C  )
% 
% [ Tab , Fi , N ] = maksttc (  w  ,  maxdt  ,  ...  )
% [ sttc , dt ] = maksttc (  Tab  ,  Fi  ,  N  ,  A  ,  B  )
% 
% MET Analysis Kit. Computes spike time tiling coefficient of Cutts and
% Eglen ( 2014 ). This provides a trial-by-trial estimation of spike train
% correlation. It is a bounded, symmetrical estimate that can distinguish
% no correlation from anti-correlation, and is robust to changes in firing
% rate and the amount of data. STTC is evaluated at all possible delta-t
% values. That is, it is evaluated at all possible time scales within the
% range of the analysis window.
% 
% w is a two-element vector defining the analysis window, where w( 1 ) is
% the starting time and w( 2 ) is the ending time of the window. STTC is
% evaluated at each possible delta-t from 0.001 to w( 2 ) - w( 1 ) in
% millisecond steps between spikes that occur within the time interval from
% w( 1 ) to w( 2 ). Hence, there are ( w( 2 ) - w( 1 ) ) / 1000 + 1
% possible delta-t values, rounded up and including 0. maxdt sets the
% maximum delta-t that may be used , in milliseconds ; if empty then all
% delta-t values are used.
%
% If there are only two neurones or spike clusters of interest, then their
% sets of spike trains can be provided in A and B. Each will contain the
% spikes from one neurone/spike-cluster or the other. A and B must be cell
% array vectors of the same length, where A{ i } and B{ i } contain the
% spike trains from the ith trial. STTC is then computed for each trial.
% Alternatively, C can be given. It must be a cell array where neurones/
% spike-clusters are indexed across columns, and trials are indexed across
% rows. For all A, B and C, each element must contain a vector of single or
% double floating point numbers that provide spike times ( in seconds ) in
% chronological order. Empty [ ] place holders can be used when there is no
% data available ; but the corresponding STTC value is undefined, so NaN
% will be returned.
% 
% sttc contains the STTC values. If A and B are given then sttc is a W x T
% matrix of STTC values. The W Delta-t values are indexed over rows ( order
% ascending ) and trials are indexed over columns. If, C is given then sttc
% becomes a W x ( M ^ 2 - M ) / 2 x T matrix, where M is the number of
% columns in C. Delta-t values are still indexed over rows. But columns are
% indexed over unique spike-cluster pairs and dimension 3 is indexed over
% trials. Cluster pairs are assigned to columns of sttc like this:
% 
%   p = 0 ;
%   for  a = 1 : M - 1
%     for  b = a + 1 : end
%       p = p + 1 ;
%       sttc( : , p , : ) = pair p's STTC at all delta-t on all trials
%     end
%   end
% 
% sttc will be NaN any time that at least one spike train was an empty
% matrix. dt will return all delta-t values in register with the rows of
% sttc. sttc and dt contain single precision floating point numbers.
% 
% If exactly three output arguments are requested then these are
% interpreted to be Tab, Fi, and N. Tab contains the proportion of time
% within delta-t of each spike. This is a single floating point matrix that
% indexes delta-t over rows, neurones/spike-clusters over columns, and
% trials over the 3rd dimension. Fi contains the index of the first spike
% in each trial that is within the analysis window defined by w. It is a
% uint32 matrix indexing neurones/spike-clusters over rows and trials over
% columns. Likewise, N is the number of spikes within the analysis window
% on each trial, and it has the same type and configuration as Fi. The
% purpose of this output is that it can be precomputed for a large number
% of neurones/spike-clusters. Then, STTC for each pair of neurones/clusters
% can be computed without re-computing Tab, Fi, and N each time. This is
% necessary for large data sets when there is insufficient memory to
% compute all pairwise correlations in a single call to maksttc.
% 
% The final form of maksttc allows precomputed Tab, Fi, and N to be given
% as input arguments. STTC is computed for single neurone/spike-cluster
% pair by accessing the subset of these matrices, and then providing A and
% B as normal. For example, if a and b are the neurone/cluster indices for
% a single pair, then STTC is found by:
% 
%   [ Tab , Fi , N ] = maksttc( w , maxdt , C ) ;
%   A = C( : , a ) ;  B = C( : , b ) ;
%   i = [ a , b ] ;
%   sttc = maksttc( Tab( : , i , : ) , Fi( i , : ) , N( i , : ) , A , B ) ;
% 
% Given Tab, Fi, and N, w is uneccesary and maxdt is inferred. As input
% arguments, Tab must have 2 columns while Fi and N must have 2 rows, each.
% 
% 
% Algorithm:
% 
% An O( n ) algorithm is used here. In other words, the number of
% computations required is at most a linear function of the size of the
% input. A simple implementation of STTC could be O( n ^ 3 ) i.e. the
% number of computations would be at most a cubic function of the size of
% the input. The O( n ) algorithm makes the assumption that spike times are
% in chronological order ; this is a minimal assumption because spikes
% occur in chronological order and are typically recorded in the order that
% they occurred. STTC is computed by:
% 
%  STTC = 0.5  *  ( ( Pa - Tb )/( 1 - PaTb )  +  ( Pb - Ta )/( 1 - PbTa ) )
% 
% For a given value of delta-t. If a and b are each spike trains with Na
% and Nb spikes times each in seconds, then Ta is the proportion of time
% within delta-t seconds of any spike in a, and Pa is the proportion of
% spikes from a that are within delta-t seconds of any spike in b ; Tb and
% Pb are the same measures again for spike train b.
%
% Let w define a time window in which w( 1 ) <= min( min( a ) , min( b ) )
% and w( 2 ) >= max( max( a ) , max( b ) ). ceil( ( w( 2 ) - w( 1 ) ) /
% 1000) is the number of W delta-t values at millisecond time steps, where
% the function ceil( x ) rounds x up to the next integer value. Hence,
% delta-t values range from 0 to ( W - 1 ) / 1000 seconds.
% 
% To calculate Ta ( and Tb ) , we start by summing the delta-t seconds on
% either side of each spike in a. Thus:
% 
%   Ta( dt ) = 2 * dt * Na / W , where dt is a value from 1 to W
% 
% But this is valid only so long as dt / 1000 is less than half of the
% shortest inter-spike-interval (ISI) between spikes in a. Once it equals
% and surpasses the shortest ISI, we compensate by realising that one ISI
% involves two spikes. Hence we must subtract 2 * dt, then add back the
% minimum ISI. Generically:
% 
%   Ta( dt ) = (  2( Na - K ) * dt  +  sum(k = 1 to K) isi(k)  )  /  W
% 
% where isi is the set of Na - 1 interspike intervals from a in
% milliseconds, sorted acsending. K is the number of ISI values less than
% 2 * dt. sum( i = 1 : N ) is used here as a substitute for sigma notation,
% indicating a sum of N values indexed starting from 1. The final
% consideration is when dt is greater than or equal to the time between
% w( 1 ) and a( 1 ) or between a( Na ) and w( 2 ). Analogous to ISI's, we
% must subtract one dt for each end that is covered, then add the duration
% back in:
% 
%   Ta( dt ) = ( N * dt  +  s( dt )Ts  +  e( dt )Te  +  Sisi )  /  W
%   
% where N = 2( Na - K ) - s( dt ) - e( dt ) and Sisi = sum(k = 1 to K)
% isi(k). Ts = a( 1 ) - w( 1 ) and Te = w( 2 ) - a( Na ). s( dt ) = 1 if
% Ts < dt, and s( dt ) = 0 otherwise. Similarly, e( dt ) = 1 if Te < dt,
% and zero otherwise.
% 
% To compute Ta in O( n ), we can use each ISI value to determine what dt
% value will surpass it. Hence, when dt >= isi(k) / 2 then we add 1 more to
% K. The first such dt' = ceil( isi(k) / 2 ). If we make K and Sisi each a
% set of W values initialised to zero, then we can use the dt' to
% accumulate the number of intervals surpassed by any given dt value:
% 
%   for k = 1 to Na - 1
%     dt = ceil( isi(k) / 2 )
%     K( dt ) = K( dt )  +  1
%     Sisi( dt ) = Sisi( dt )  +  isi(k)
%   end
% 
% The cumulative sum of K and Sisifrom 1 to W will then provide the number
% of ISIs surpassed by any given dt, while the cumulative sum of Sisi will
% give the equivalent value of sum(k = 1 to K(dt)) isi(k):
% 
%   for i = 2 to W
%     K( i ) = K( i )  +  K( i - 1 )
%     Sisi( i ) = Sisi( i )  +  Sisi( i - 1 )
%   end
%  
% Thus, each Ta becomes:
% 
%   for i = 1 to W
%     Ta( i ) = ( N( i ) * dt +  s( dt )Ts +  e( dt )Te  +  Sisi( i ) ) / W
%   end
%   
%  where N( i ) = 2( Na - K( i ) ) - s( i ) - e( i ). Tb is computed the
%  same way from b.
% 
% Pa and Pb are easier to compute. This will require computing indices with
% the function:
% 
%   F( x , y ) = ceil(  ( x - y )  /  0.001  )
% 
% Use index variables ia and ib to track the current spike in a and b. Both
% start at 1. We begin by accumulating data from the head of the spike
% train that starts first:
% 
%   if  a( 1 ) <= b( 1 )
% 
%     while  a( ia ) <= b( 1 )
%       dt = F( b( 1 ) , a( ia ) ) 
%       Pa( dt ) = Pa( dt ) + 1
%       ia = ia + 1
%     end
% 
%   else
% 
%     while  b( ib ) <= a( 1 )
%       dt = F( a( 1 ) , b( ib ) ) 
%       Pb( dt ) = Pb( dt ) + 1
%       ib = ib + 1
%     end
% 
%   end
% 
% At this point, the next chronological spike is b( 1 ) in the first case
% and a( 1 ) in the second. It is guaranteed to have at least one spike
% from the other train before it. This allows us to accumulate data from
% spikes with neighbours from the other train on both sides:
% 
%   while  ia <= Na AND ib <= Nb
% 
%     if  a( ia ) <= b( ib )
% 
%       d1 = F( a( ia ) , b( ib - 1 ) )
%       d2 = F( b( ib ) , a( ia ) )
%       dt = min( d1 , d2 )
%       Pa( dt ) = Pa( dt ) + 1
%       ia = ia + 1
% 
%     else
% 
%       d1 = F( b( ib ) , a( ia - 1 ) )
%       d2 = F( a( ia ) , b( ib ) )
%       dt = min( d1 , d2 )
%       Pb( dt ) = Pb( dt ) + 1
%       ib = ib + 1
% 
%     end
% 
%   end while
% 
% At last, we accumulate data from tailing spikes in one train that all
% come after the final spike in the other. We can tell which is the lagging
% train using ia and ib following the break condition for the previous
% step:
% 
%   if  ia <= Na
%
%     while  ia <= Na
%       dt = F( a( ia ) , b( Nb ) )
%       Pa( dt ) = Pa( dt )  +  1
%       ia = ia  +  1
%     end
%
%   elseif  b <= Nb
%
%     while  ib <= Nb
%       dt = F( b( ib ) , a( Na ) )
%       Pb( dt ) = Pb( dt )  +  1
%       ib = ib  +  1
%     end
%
%   end
% 
% Pa( dt ) now give the number of spikes from a that are dt milliseconds
% from a neighbouring spike in b. Pa( dt ) has the equivalent for b. Again,
% cumulative sums are used to make Pa( dt ) equal the number of spikes in a
% that are within dt milliseconds of a spike in b. To further normalise by
% Na gives us the proportions we need to compute STTC:
% 
%   for i = 2 to W
%     Pa( i ) = ( Pa( i ) + Pa( i - 1 ) )
%     Pb( i ) = ( Pb( i ) + Pb( i - 1 ) )
%     Pa( i - 1 ) = Pa( i - 1 ) / Na
%     Pb( i - 1 ) = Pb( i - 1 ) / Nb
%   end
%   
%   Pa( W ) = Pa( W ) / Na
%   Pb( W ) = Pb( W ) / Nb
% 
% At last, STTC is computed at each delta-t by:
% 
%   for dt = 1 to W
%     STTC( dt ) = 0.5  *  (  ( Pa(dt) - Tb(dt) )/( 1 - Pa(dt)*Tb(dt) )  +
%       ( Pb(dt) - Ta(dt) )/( 1 - Pb(dt)*Ta(dt) )  ) 
%   end
% 
% 
% Implementation:
% 
% Uses Matlab's Parallel Computing Toolbox.
% 
% 
% Reference:
% 
% Cutts CS, Eglen SJ. 2014. Detecting Pairwise Correlations in Spike
%   Trains: An Objective Comparison of Methods and Application to the Study
%   of Retinal Waves. J Neurosc, 34(43):14288-14303.
% 
% 
% Written by Jackson Smith - March 2018 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % Number of inputs for ( ... , A , B ) function call format
  NARGAB = 4 ;
  
  % Number of inputs for ( Tab , Fi , N , A , B ) function call format
  NINTFN = 5 ;
  
  % Minimum time step in seconds , one millisecond
  MINSTP = 0.001 ;
  
  
  %%% Check input %%%
  
  % Check number of input arguments
   narginchk (  3  ,  NINTFN  )
  nargoutchk (  0  ,  3  )
  
  % The number of inputs tells us what form we have. If less than 5 then we
  % have w and maxdt.
  if  nargin  <  NINTFN
    
    % Name inputs
    [ w , maxdt ] = varargin{ 1 : 2 } ;
  
    % Duration of analysis window
    W = round(  diff(  w  )  ,  6  ) ;

    % Check w
    if  numel( w )  ~=  2  ||  ~ isnumeric( w )  ||  ~ isreal( w )  ||...
          w( 2 )  <=  w( 1 )  ||  W  <  MINSTP

      error (  'MAK:maksttc:w'  ,  [ 'maksttc: w must be 2-element ' , ...
        'real value vector where w( 1 ) < w( 2 ) and %f <= ' , ...
          'w( 2 ) - w( 1 )' ]  ,  MINSTP  )

    end % check w

    % Finish computing number of milliseconds spanned by w
    W = ceil (  W  *  1000  ) ;

    % Check maxdt , if empty then it is ignored
    if  ~ isempty (  maxdt  )

      % Must be a scalar real number of at least the minimum delta-t and up
      % to the maximum
      if  ~ isscalar (  maxdt  )  ||  ~ isnumeric (  maxdt  )  ||  ...
          ~ isreal (  maxdt  )  ||  maxdt < 1e3 * MINSTP  ||  W < maxdt

        error (  'MAK:maksttc:maxdt'  ,  [ 'maksttc: maxdt must be ' , ...
          'a scalar real number in the range of %d to %d' ]  ,  ...
            1e3 * MINSTP  ,  W  )

      end

      % Take minimum value for W
      W = min (  maxdt  ,  W  ) ;

    end % check maxdt

    % Add 1 to W to include delta-t of 0 and W becomes the number of delta-t
    % values
    W = W + 1 ;
  
  end % less than 5 input arg forms
  
  % A and B form of input
  if  nargin  >=  NARGAB
    
    % Get A and B
    [ A , B ] = varargin{ end - 1 : end } ;
    
    % Make sure that these are cell arrays
    if  ~ iscell (  A  )  ||  ~ iscell (  B  )
      
      error (  'MAK:maksttc:ABcell'  ,  ...
        'maksttc: A and B must be cell arrays'  )
      
    % Must be vectors
    elseif  ~ isvector (  A  )  ||  ~ isvector (  B  )
      
      error (  'MAK:maksttc:ABvect'  ,  ...
        'maksttc: A and B must be vectors'  )
      
    % Must have equal number of elements
    elseif  numel (  A  )  ~=  numel (  B  )
      
      error (  'MAK:maksttc:ABnum'  ,  ...
        'maksttc: A and B must have the same number of elements'  )
      
    end % check A and B
    
    % Build C
    C = [  A( : )  ,  B( : )  ] ;
    
  % C form of input
  else
    
    % Get C
    C = varargin{ end } ;
    
    % Must be a cell array
    if  ~ iscell (  C  )
      
      error (  'MAK:maksttc:Ccell'  ,  'maksttc: C must be a cell array'  )
      
    % Must not exceed 2 dimensions and must not be empty
    elseif  ~ ismatrix (  C  )  ||  isempty (  C  )
      
      error (  'MAK:maksttc:Ccell'  ,  ...
        'maksttc: C must not be empty or exceed 2 dimensions'  )
      
    % Must not have fewer than 2 columns
    elseif  size (  C   ,  2  )  <  2
      
      error (  'MAK:maksttc:Cncols'  ,  ...
        'maksttc: C not have fewer than 2 columns'  )
      
    end % check C
    
  end % get input
  
  % Check elements of C
  parfor  i = 1 : numel( C )
    
    % This spike train
    c = C{ i } ;
    
    % Check that it has all these attributes
    X( i ) = ( isvector( c ) || isempty( c ) )  &&  ...
    isnumeric( c )  &&  ~isinteger( c )  &&  isreal( c ) ;
  
  end
          
  % All spike trains must have listed attributes
	if  any (  ~ X( : )  )
    
    error (  'MAK:maksttc:spktrains'  ,  [ 'maksttc: all spike ' , ...
      'trains must be real-numbered vectors of type single or ' , ...
        'double , or empty' ]  )
      
  else
    
    clear  X
    
  end % check elements of C
  
  % Number of trials and spike clusters
  [ Nt , Ns ] = size (  C  ) ;
  
  % 5 input argument form
  if  nargin  ==  NINTFN
    
    % There must be no more than 2 output arguments
    if  2  <  nargout
      
      error (  'MAK:maksttc:argsout'  ,  [ 'maksttc: too many ' , ...
        'outputs requested from %d input argument form of function ' , ...
          'call' ]  ,  NINTFN  )
      
    end
    
    % Check data
    X = {  {  'T'  ,  'single'  ,  3  }  ;
           {  'Fi' ,  'uint32'  ,  2  }  ;
           {  'N'  ,  'uint32'  ,  2  }  } ;
         
    % Input args
    for  i = 1 : numel (  X  )
      
      % Name and type
      [ n , t , reqnd ] = X{ i }{ : } ;
      
      % Number of dimensions
      nd = ndims (  varargin{ i }  ) ;
      
      % Size of argument
      s = size (  varargin{ i }  ) ;
      
      % Check type
      if  ~ isa (  varargin{ i }  ,  t  )
        
        error (  'MAK:maksttc:argtype'  ,  ...
          'maksttc: %s must be of type %s'  ,  n  ,  t  )
        
      % Check number of dimensions
      elseif  nd  ~=  reqnd
        
        error (  'MAK:maksttc:argdim'  ,  ...
          'maksttc: %s must have %d dimensions'  ,  n  ,  reqnd  )
        
      % Check number of neurones/spike-clusters
      elseif  s (  nd - 1  )  ~=  2
        
        error (  'MAK:maksttc:argnumclusts'  ,  ...
          'maksttc: %s must have length 2 in dimension %d'  ,  ...
            n  ,  nd - 1  )
          
      % Check number of trials
      elseif  s (  end  )  ~=  Nt
        
        error (  'MAK:maksttc:argnumtrials'  ,  ...
          'maksttc: %s must have %d trials'  ,  n  ,  Nt  )
        
      end % check arg
      
    end % args
    
    % Give names to input args
    [ T , Fi , N ]  = varargin{ 1 : 3 } ;
    
    % Infer number of delta-t values
    W = size (  T  ,  1  ) ;
    
  end % 5 arg form
  
  
  %%% Preparation %%%
  
  % It is actually more convenient if C is arranged as spike-clusters over
  % rows by trials over columns. Transpose.
  C = C' ;
  
  % Fi and N do not exist if fewer than 5 input arguments were given
  if  nargin  <  NINTFN
    
    % Enumerate each millisecond time-scale
    dt = enumdt (  W  ) ;

    % Find all spikes within analysis window , return first and last spike
    % indices for each cluster and trial
    parfor  i = 1 : numel (  C  )

      [ Fi( i ) , N( i ) ] = getspks (  C{ i }  ,  w  ) ;

    end

    % Reshape to match C
    Fi  = reshape (  Fi   ,  size(  C  )  ) ;
    N   = reshape (  N    ,  size(  C  )  ) ;
    
    % Compute proportion of time for all spike trains
    parfor  i = 1 : numel (  C  )
    
      T( : , i ) = getT (  w ,  W ,  dt ,  C{ i } ,  Fi( i ) ,  N( i )  ) ;

    end

    % T is now a delta-t by spike-clusters x trials 2D matrix. Make it a
    % delta-t by spike-clusters by trials 3D matrix.
    T = reshape (  T  ,  [  W  ,  Ns  ,  Nt  ]  ) ;
    
    % Return Tab, Fi, and N
    if  nargout  ==  3

      varargout( 1 : 3 ) = {  T  ,  Fi  ,  N  } ;
      return

    end
  
  end % make Fi and N
  
  
  %%% Indexing spike-cluster pairs %%%
  
  % Figure out the number of unique cross-correlated spike cluster pairs
  Np = ( Ns ^ 2  -  Ns )  /  2 ;
  
  % Determine the linear index of each pair in the lower-triangular part of
  % a Ns x Ns matrix
  I = find (  tril (  true (  Ns  )  ,  -1  )  )' ;
  
  % Translate this into cluster indices i.e. columns of C. But, from the
  % perspective of the Ns x Ns matrix, these are row and column sub-
  % scripts.
  [ bi , ai ] = ind2sub (  [ Ns , Ns ]  ,  I  ) ;
  
  
  %%% Compute STTC %%%
  
  % Implicit pointers to T. Since neither Ta nor Tb will be assigned to,
  % copy-on-write will make them both reference T when accessed.
  Ta = T ;
  Tb = T ;
  
  % Allocate output variable. parfor output gets placed into this
  sttc = zeros (  W  ,  Np  ,  Nt  ,  'single'  ) ;
  
  % Pairs of spike clusters
  for  p = 1 : Np
    
    % Row and column of a Ns x Ns matrix , or spike cluster indices of a
    % unique pair
    a = ai( p ) ;  b = bi( p ) ;
    
    % Access data for cluster A
    A =  C( a , : ) ;
  Fia = Fi( a , : ) ;
   Na =  N( a , : ) ;
    
    % Access data for cluster B
    B =  C( b , : ) ;
  Fib = Fi( b , : ) ;
   Nb =  N( b , : ) ;
      
    % Trials
    parfor  i = 1 : Nt

      % Proportion of spikes within delta-t of each other
      [ Pa , Pb ] = getP (  W ,  A{ i } ,  B{ i } ,  Fia( i ) ,  ...
        Fib( i ) ,  Na( i ) ,  Nb( i )  ) ;

      % Compute STTC for this trial. Due to the finite precision of
      % computers, we must use the min function to get rid of NaN values
      % that occur due to division by zero. This results when delta-t
      % grows large and P* == T == 1.
      sttc( : , p , i ) = 0.5  *  ( ...
   min(  ( Pa  -  Tb(:,b,i) )  ./  ( 1  -  Pa .* Tb(:,b,i) ) ,  1  )  + ...
   min(  ( Pb  -  Ta(:,a,i) )  ./  ( 1  -  Pb .* Ta(:,a,i) ) ,  1  )  ) ;

    end % trials
    
  end % pairs
  
  
  %%% Return arguments %%%
  
  % A B input arguments require that sttc is delta-t over rows and trials
  % over columns
  if  nargin  >=  NARGAB  ,  sttc = reshape (  sttc  ,  W  ,  Nt  ) ;  end
  
  % Return sttc
  varargout{ 1 } = sttc ;
  
  % If dt is requested then return a single
  if  1  <  nargout
    
    % dt not created yet
    if  ~ exist (  'dt'  ,  'var'  )  ,  dt = enumdt (  W  ) ;  end
    
    % Return dt
    varargout{ 2 } = single (  dt  ) ;
  
  end % ret dt
  
  
end % maksttc


%%% Subroutines %%%

% Enumerate delta-t values
function  dt = enumdt (  W  )
  
  dt = (  0 : W - 1  )' ;
  
end % enumdt


% Index of first and last spike within specified time window
function  [ nf , n ] = getspks (  s  ,  w  )
  
  % Locate all spikes in window
  i = w( 1 )  <=  s  &  s  <=  w( 2 ) ;
  
  % Count spikes in window
  n = uint32 (  sum(  i  )  ) ;
  
  % Spikes found in window , find index of first spike in window
  if  n
    
    nf = uint32 (  find(  i  ,  1  ,  'first'  )  ) ;
  
  % No spikes in window , return placeholder value
  else
    
    nf = zeros (  1  ,  'uint32'  ) ;
    
  end
  
end % getspks


% Computes delta-t in milliseconds between two times. Adds 1 to account for
% delta-t of zero.
function  d = F (  x  ,  y  )

  d = ceil (  (  x  -  y  )  *  1000  )  +  1 ;
  
end % F

% Compute proportion of time within any of N spikes from train A for each
% millisecond time scale from 1 to W
function  T = getT (  w  ,  W  ,  dt  ,  A  ,  Nf  ,  Na  )
  
  % Use double precision floating point for calculations
  if  ~ isa (  A  ,  'double'  )  ,  A = double (  A  ) ;  end
  
  % Spike train index vector
  ia = Nf : Nf + Na - 1 ;

  % No spikes. STTC is undefined so return NaN.
  if  Na  ==  0
    T = nan (  W  ,  1  ,  'single'  ) ;
    return
  end
  
  % Allocate number of surpassed inter-spike-intervals and accumulation of
  % surpassed ISI's. Also count the number of window start/end to spike
  % intervals surpassed by each delta-t.
  K = zeros (  W  ,  1  ) ;
  Sisi = zeros (  W  ,  1  ) ;
  
  % Inter-spike-intervals (rounded to nearest nanosecond) in milliseconds
  isi = diff(  A( ia )  )  *  1000 ;
  
  % Calculate delta-t that surpasses each isi. Since this is index for K
  % and S, we call it i. We must add 1 to account for delta-t of zero.
  i = ceil (  isi  /  2  )  +  1 ;
  
  % Accumulate isi data
  for  j = 1 : Na - 1
    
    % Get delta-t that surpasses the jth ISI
    d = i( j ) ;
    
    % Greater than the maximum allowable delta-t
    if  W  <  d  ,  continue  ,  end
    
    % Count number of ISI's surpassed at this delta-t , and accumulate
    % their sum
    K( d ) = K( d )  +  1 ;
    Sisi( d ) = Sisi( d )  +  isi( j ) ;
    
  end % surpassed ISI
  
  % Cumulative sums render K and S* into number of surpassed intervals ,
  % and sum of all surpassed ISI's at each delta-t
  K = cumsum (  K  ) ;
  Sisi = cumsum (  Sisi  ) ;
  
  % Time from start of window to first spike and from the last spike to end
  % of the window , in milliseconds
  Ts = (  A( Nf )  -  w(  1  )        )  *  1000 ;
  Te = (  w( 2 )   -  A( ia( end ) )  )  *  1000 ;
  
  % Determine which delta-t's exceed the time to first spike or time from
  % last spike to the start and end of the analysis window
  s = Ts  <  dt ;
  e = Te  <  dt ;
  
  % Compute number of delta-t intervals to sum
  N = 2 * (  double( Na )  -  K  )  -  s  -  e ;
  
  % Compute proportion of time covered at each delta-t
  T = (  N .* dt  +  s * Ts  +  e * Te  +  Sisi  )  /  ...
    (  diff( w )  *  1000  ) ;
  
  % Return single-precision floating point numbers
  T = single (  T  ) ;
  
end % getT


% Compute proportion of spikes in A within each delta-t of spikes from B ,
% and vice-versa
function  [ Pa , Pb ] = getP (  W  ,  A  ,  B  ,  a  ,  b  ,  Na  ,  Nb  )
  
  % If either train is empty then STTC is undefined so return NaN
  if  Na  ==  0  ||  Nb  ==  0
    Pa = nan (  W  ,  1  ,  'single'  ) ;
    Pb = nan (  W  ,  1  ,  'single'  ) ;
    return
  end
  
  % Index of the last spike in each train within window
  La = a + Na - 1 ;
  Lb = b + Nb - 1 ;
  
  % Allocate vectors for output , initialised to zero
  Pa = zeros (  W  ,  1  ) ;
  Pb = zeros (  W  ,  1  ) ;
  
  % Leading spikes in train A
  if  A( a )  <=  B( b )
    
    % Count spikes in A that are an exact distance to the first spike in B
    while  a  <=  La  &&  A( a )  <=  B( b )
      
      % Find delta-t separating these two spikes
      d = F (  B( b )  ,  A( a )  ) ;
      
      % Increment to next spike in A
      a = a  +  1 ;
      
      % Bigger than maximum delta-t , carry on
      if  W < d  ,  continue  ,  end
      
      % Count one more spike in A at this distance from a neighbour in B
      Pa( d ) = Pa( d )  +  1 ;
      
    end % count A near B
    
  % Leading spikes in train B
  else
    
    while  b  <=  Lb  &&  B( b )  <=  A( a )
      d = F (  A( a )  ,  B( b )  ) ;
      b = b  +  1 ;
      if  W < d  ,  continue  ,  end
      Pb( d ) = Pb( d )  +  1 ;
    end
    
  end % leading spikes
  
  % Spikes with neighbours from other train on both sides
  while  a  <=  La  &&  b  <=  Lb
    
    % Next chronological spike is from train A
    if  A( a )  <=  B( b )
      
      % Get the number of milliseconds from this spike in A to the two
      % neighbouring spikes in B that come before and after the spike in A
      d1 = F (  A( a )  ,  B( b - 1 )  ) ;
      d2 = F (  B( b )  ,  A( a )      ) ;
      
      % Get the minimum distance to a neighbouring spike
      d = min (  d1  ,  d2  ) ;
      
      % Go to next spike in A
      a = a  +  1 ;
      
      % Bigger than maximum delta-t , carry on
      if  W < d  ,  continue  ,  end
      
      % Count one more spike from A that is d milliseconds from a
      % neighbouring spike in B
      Pa( d ) = Pa( d )  +  1 ;
      
    % Next chronological spike is from train B
    else
      
      d1 = F (  B( b )  ,  A( a - 1 )  ) ;
      d2 = F (  A( a )  ,  B( b )      ) ;
      d = min (  d1  ,  d2  ) ;
      b = b  +  1 ;
      if  W < d  ,  continue  ,  end
      Pb( d ) = Pb( d )  +  1 ;
      
    end % next chronological spike
    
  end % spikes with neighbours
  
  % Trailing spikes in train A
  if  a  <=  La
    
    % Count spikes in A that are an exact distance to the last spike in B
    while  a  <=  La
      
      % delta-t that separates spike pair
      d = F (  A( a )  ,  B( Lb )  ) ;
      
      % Increment to next spike in A
      a = a  +  1 ;
      
      % Bigger than maximum delta-t , carry on
      if  W < d  ,  continue  ,  end
      
      % Count one more spike in A that is d milliseconds from a spike in B
      Pa( d ) = Pa( d )  +  1 ;
      
    end % count A near B
    
  % Trailing spikes in B
  elseif  b  <=  Lb
    
    while  b  <=  Lb
      d = F (  B( b )  ,  A( La )  ) ;
      b = b  +  1 ;
      if  W < d  ,  continue  ,  end
      Pb( d ) = Pb( d )  +  1 ;
    end
    
  end % trailing spikes
  
  % Cumulative sums give the number of spikes in one train within each
  % delta-t of a spike in the other train. Normalising by number of spikes
  % per train gives the proportions.
  Pa = single (  cumsum (  Pa  )  /  double( Na )  ) ;
  Pb = single (  cumsum (  Pb  )  /  double( Nb )  ) ;
  
end % getP

