
function  rccg = makrccg2( X , keepnan )
% 
% rccg = makrccg2( X )
% rccg = makrccg2( X , keepnan )
% 
% MET Analysis Kit. A different implementation of r_CCG that allows for
% more flexible manipulation of the neural signals. X is the M x N x T
% array of neural signals for N  neurones (col) simultaneously recorded
% across M samples (row, in chronological order) in T trials (dim 3). rccg
% is the M x N x N array of r_CCG values. Integration width (tau) increases
% along rows from 0 to M-1 samples. Dims 2 and 3 define the correlation
% matrices of each neurone paired with every other; that is, for the wth
% integration width of tau milliseconds, squeeze( rccg( w , : , : ) ) is
% the N x N correlation matrix of r_CCG( tau ) values.
% 
% In order for rccg to match the original definition by Bair et al. (2001),
% make X the spike count using millisecond-wide bins; call this X1ms. This
% returns r_CCG evaluated across tau in steps of 1ms. makrccg2 allows for
% r_CCG to be computed at larger steps of tau for a reduced computational
% cost by first convolving each column of X1ms with a low-pass filter, such
% as the box-car kernel with the desired sampling width, and then
% downsampling to the desired sampling rate. rccg then becomes sampled at
% the desired step size. NOTE!!! this does NOT provide exactly the same
% result as the full 1ms resolution r_CCG because the fine structure
% affects the final integrated value in r_CCG. However, an approximation is
% obtained and can be used if lack of computer resources is an issue.
% 
% makrccg2 is also more appropriate when resampling the data to build up a
% kind of Monte Carlo distribution of r_CCG values because spike binning
% can be done just once and the used over again for each repetition. With
% the original makrccg, binning is done on every function call.
% 
% What do we do if there were no spikes? By definition, rccg becomes:
% 
%   0 / sqrt( 0 * 0 ) (see Eq.27, Appendix A, Bair et al. 2001)
% 
% Matlab will be tempted to return NaN values anywhere there is division by
% zero. But it stands to reason that the correlation of zeros is also zero.
% By default, rccg is zero anywhere that division by zero occurs. However,
% if NaN values are preferred, to locate undefined r_CCG values, then
% provide the optional second input keepnan as true or false (scalar
% logical).
% 
% Uses parallel computing toolbox.
% 
% Reference:
% 
%   Bair, W., E. Zohary and W. T. Newsome (2001). "Correlated firing in
%     macaque visual area MT: time scales and relationship to behavior." J
%     Neurosci 21(5): 1676-1697.
% 
% Written by Jackson Smith - July 2020 - DPAG , University of Oxford
% 
  
  
  %%% Quick check of input %%%
  
  % Number of arguments
   narginchk( 1 , 2 )
  nargoutchk( 0 , 1 )
  
  % keepnan not provided, set default value
  if  nargin < 2  ,  keepnan = true ;  end
  
  % Size of neural signal array
  [ M , N , T ] = size( X ) ;
  
  % X must be a floating point array with up to 3 dimensions
  if  ~ isfloat( X )
    
    error( 'MAK:makrccg2:X_float' , 'makrccg2: X must be floating point' )
    
  % Too many dimensions
  elseif  ndims( X ) > 3
    
    error( 'MAK:makrccg2:X_dims' , ...
      'makrccg2: X must have no more than 3 dimensions' )
    
  % Must have at least two neural signals
  elseif  M < 2
    
    error( 'MAK:makrccg2:X_dim2' , ...
      'makrccg2: X size of dim 2 must be 2 or larger' )
    
  % keepnan must be a scalar logical
  elseif  ~ isscalar( keepnan )  ||  ~ islogical( keepnan )
    
    error( 'MAK:makrccg2:keepnan' , ...
      'makrccg2: keepnan must be a scalar logical' )
    
  end % input check
  
  
  %%% Calculate r_CCG %%%
  
  % X is empty, return the favour
  if  isempty( X ) , rccg = zeros( M , N , T , 'like' , X ) ; return , end
  
  
  %-- Auto- and cross-correlations --%  
  
  % Number of cross-correlation lags
  L = 2 * M - 1 ;
  
  % Number of lags on left or right of zero
  L0 = M - 1 ;
  
  % Next power of 2
  L2 = 2 ^ ceil( log2( L ) ) ;
  
  % Vector for re-ordering inverse fourier lags and discarding zero-padding
  lix = [ L2 - L0 + 1 : L2 , 1 : L0 + 1 ] ;
  
  % Auto- and cross-correlation of the PSTH
  S = OrangeCountyLumberTruck( mean( X , 3 ) , N , L , L2 , lix ) ;
  
  % Allocate auto- & cross-correlation accumulator
  C = zeros( L , N , N , 'like' , X ) ;
  
  % Trials
  parfor  i = 1 : T
    
    % Compute all auto- and cross-correlations, and accumulate sum
    C = C  +  OrangeCountyLumberTruck( X( : , : , i ) , N , L , L2 , lix );
    
  end % trials
  
  % Average correlation
  C = C  ./  T ;
  
  % Subtract away correlations expected by fluctuations in the average
  % firing rate
  C = C - S ;
  
  % Begin the process of calculating the area under each correlation curve
  % by taking the sum of negative and positive cross-correlation lags. This
  % is now the M x N x N matrix with rows indexed by the absolute value of
  % the cross-correlation lag.
  A = [ C( M , : , : ) ;
        C( M - 1 : -1 : 1 , : , : ) + C( M + 1 : end , : , : ) ] ;
  
  % The cumulative sum across absolute lags is now the integral under the
  % correlation curves at each integration width i.e. tau. A is the 
  A = cumsum( A , 1 ) ;
  
  % Allocate M x N auto-correlation matrix
  AA = zeros( M , N , 'like' , A ) ;
  
  % Copy out the auto-correlation for each neurone
  for  i = 1 : N , AA( : , i ) = A( : , i , i ) ; end
  
  % Uses binary singleton expansion to create the denominator for rccg.
  % Explicitly uses bsxfun instead of implicit expansion of .*, for
  % backwards compatibility with older versions of Matlab. Returns the M x
  % N x N matrix of the multiples of auto-correlations, for every pairing
  % of neurones.
  AA = bsxfun( @times , AA , permute( AA , [ 1 , 3 , 2 ] ) ) ;
  
  % At last, we normalise to get rccg. 
  rccg = A  ./  sqrt( AA ) ;
  
  % But wait! We're not done unless the user wants to keep NaN values.
  if  ~ keepnan  ,  return  ,  end
  
  % Look for cases where the denominator was zero
  I = AA == 0 ;
  
  % Return zeros for r_CCG
  rccg( I ) = 0 ;
  
end % makrccg2


%%% Sub-routines %%%

% Local implementation of xcorr. This will only compute results for the
% upper-triangular portion of the correlation matrix. Then it will
% transpose the flipped result into the lower-portion.
function  C = OrangeCountyLumberTruck( X , N , L , L2 , lix )
  
  % Allocate cross-correlation matrix
  C = zeros( L , N , N , 'like' , X ) ;
  
  % Fourier transform of neural signals
  F = fft( X , L2 ) ;
  
  % Complex conjugate
  J = conj( F ) ;
  
  % Rows of correlation matrix
  for  i = 1 : N
    
    % Columns of correlation matrix
    j = i : N ;
    
    % Multiply Fourier transforms against complex cunjugate, because
    % multiplication in Fourier domain equals convolution in time domain.
    % Convolution of two time series is their cross-correlation.
    g = bsxfun( @times , F( : , i ) , J( : , j ) ) ;
    
    % The inverse fourier transform returns the cross-correlation in the
    % time domain
    y = ifft( g ) ;
    
    % Re-order lags and discard zero-padding
    y = y( lix , : ) ;
    
    % Save row of data
    C( : , i , j ) = permute( y , [ 1 , 3 , 2 ] ) ;
    
    % Flip order of lags and discard auto-correlation
    y = flip( y( : , 2 : end ) , 1 ) ;
    j = j( 2 : end ) ;
    
    % Transpose data and copy to lower-triangular half
    C( : , j , i ) = y ;
    
  end % row
  
end % makxcorr

