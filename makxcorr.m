
function  [ y , lags ] = makxcorr (  a  ,  b  )
% 
% y = makxcorr (  a  )
% y = makxcorr (  a  ,  b  )
% [ y , lags ] = makxcorr ( ... )
% 
% MET Analysis Kit. Computes cross-correlation. Matlab's xcorr can run
% slowly for large inputs because it uses (at least, as of R2015b) for
% loops. makxcorr makes an increased use of vectorisation in order to run
% faster; otherwise, it is based on the same principal in which signals in
% the time domain are converted to the Fourier domain, where they are
% multiplied together before being converted back to the time domain.
% 
% Input arguments must be floating point matrices. It is assumed that
% samples are indexed across rows. Columns are treated as separate signals.
% 
% With one input, the matrix auto-correlation is calculated such that the
% correlation of all pairings of columns are returned in y. If a is an M by
% N matrix, then y is the 2M-1 by N by N matrix with lags indexed over
% rows. For the Lth lag, sy = squeeze( y( L , : , : ) ) is the square
% correlation matrix such that sy( i , j ) is the correlation of the input
% signals a( : , i ) and a( : , j ).
% 
% When two input arguments a and b with sizes M x N and P x Q are given
% then y becomes the M+P-1 by N by Q cross-correlation matrix. For the Lth
% lag, sy = squeeze( y( L , : , : ) ) is the N by Q matrix such that
% sy( i , j ) is the correlation of signals a( : , i ) and b( : , j ).
% 
% Written by Jackson Smith - August 2019 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Two args
  twoarg = 1  <  nargin ;
  
  % Check arg 'a'
  if  ~ isfloat (  a  )  ||  ~ ismatrix (  a  )
    
    error (  'MAK:makxcorr:a_invalid'  ,  ...
      'makxcorr: input argument ''a'' must be a floating point matrix'  )
    
  % Check arg 'b'
  elseif  twoarg  &&  (  ~isfloat (  b  )  ||  ~ismatrix (  b  )  )
    
    error (  'MAK:makxcorr:b_invalid'  ,  ...
      'makxcorr: input argument ''b'' must be a floating point matrix'  )
    
  end % check args
  
  
  %%% Prep %%%
  
  % Size of input args a ...
  [ M , N ] = size (  a  ) ;
  
  % ... and b
  if  twoarg
    
    [ P , Q ] = size (  b  ) ;
  
  else
    
    % Point to a's size for generic code
    P = M ;  Q = N ;
    
  end
  
  % Maximum number of lags , we will zero pad whichever input arg has fewer
  % rows
  nlags = 2 * max( M , P )  -  1 ;
  
  % Size of discrete fourier transform, round up to the next power of 2
  n = 2  ^  ceil ( log2(  nlags  ) ) ;
  
  
  %%% Compute cross-correlations %%%
  
  % We start by taking the fourier transform of input arg a
  f = fft (  a  ,  n  ) ;
  
  % Now get the complex conjugate, this depends on the number of inputs
  if  twoarg
    
    % b was provided, return the conjugate of its fourier transform
    c = conj ( fft(  b  ,  n  ) ) ;
    
  else
    
    % Just 'a' is given, so take conjugate of existing fourier transform
    c = conj (  f  ) ;
    
  end % complex conjugate
  
  % Reshape conjugate matrix so that columns are indexec across dim 3
  if  1  <  Q  ,  c = permute (  c  ,  [ 1 , 3 , 2 ]  ) ;  end
  
  % Multiply fourier transforms against complex cunjugates in every
  % possible combination of pairings
  g = bsxfun (  @times  ,  f  ,  c  ) ;
  
  % The inverse fourier transform returns the cross-correlation in the time
  % domain
  y = ifft (  g  ) ;
  
  % Number of lags to the left or right of lag zero
  nlags = ( nlags  -  1 )  ./  2 ;
  
  % Put lags into the correct order for return, indexing through negative
  % lags up to lag zero, and then through positive lags, in that order
  y = y( [ end - nlags + 1 : end , 1 : nlags + 1 ] , : , : )  ;
  
	% Lags requested as well
  if  1  <  nargout  ,  lags = ( -nlags : +nlags )' ;  end
  
  
end % makxcorr

