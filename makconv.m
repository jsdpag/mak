
function  c = makconv (  x  ,  k  ,  t  )
% 
% c = makconv (  x  ,  k  ,  t  )
% 
% MET Analysis Kit. Performs linear convolution on each column of matrix x
% using the kernel in vector k. If x has M rows while k has K points, then
% c is a M + K - 1 row array of convolutions such that c( : , i ) is the
% convolution of x( : , i ) and k. Convolution occurs over rows, but x can
% be multi-dimensional; c will match the size of x in dims 2 and higher.
% 
% Optional char t says what kind of kernel k is, and reshapes c
% accordingly. t can be:
% 
%   'p' - k is a predictive kernel, remove the first K - 1 rows of c
%   's' - k is a symmetric kernel, remove the first floor( K / 2 ) rows and
%     the last floor( ( K - 1 ) / 2 ) rows to return the central portion of
%     the convolution that corresponds to x.
%   'c' - k is a causal kernel, remove the last K - 1 rows of c.
% 
% If t is an empty string then it is ignored.
% 
% Written by Jackson Smith - August 2019 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Optional input not provided, set to default of empty string
  if  nargin  <  3  ,  t = '' ;  end
  
  % x is numeric
  if  ~ isnumeric (  x  )
    
    error (  'MAK:makconv:x'  ,  'makconv: x must be numeric'  )
    
  % k is numeric vector
  elseif  ~ isvector (  k  )  ||  ~ isnumeric (  k  )
    
    error (  'MAK:makconv:k'  ,  'makconv: k must be numeric vector'  )
    
  % t is char with length of 1 at most
  elseif  ~ ischar (  t  )  ||  1  <  numel (  t  )
    
    error (  'MAK:makconv:t'  ,  ...
      'makconv: t must be char with length of 1 at most'  )
    
  end % check in
  
  % Make sure that k is a column vector
  if  ~ iscolumn (  k  )  ,  k = k( : ) ;  end
  
  
  %%% Convolution %%%
  
  % Data points in each arg, M is the number of rows of x, K is number of
  % point in kernel
  M = size (  x  ,  1  ) ;
  K = numel (  k  ) ;
  
  % Number of dimensions of x, excluding rows
  N = ndims (  x  )  -  1 ;
  
  % Length of linear convolution
  L = M + K - 1 ;
  
  % Optimal size for fast fft run time
  nf = 2  ^  ceil( log2(  L  ) ) ;
  
  % Fourier transform of the inputs
  fx = fft (  x  ,  nf  ) ;
  fk = fft (  k  ,  nf  ) ;
  
  % Convolution is equal to multiplication in Fourier domain
  c = ifft ( bsxfun(  @times  ,  fx  ,  fk  ) ) ;
  
  % Get colon operator for each dimension of c, excluding rows
  co = repmat (  { ':' }  ,  1  ,  N  ) ;
  
  % Get rid of excess data points used for Fourier transformations
  c( L + 1 : end , co{ : } ) = [ ] ;
  
  % Kernel type
  switch  t
    
    % Predictive, remove first K - 1 rows
    case  'p'  ,  c( 1 : K - 1 , co{ : } ) = [ ] ;
      
    % Symmetric, remove the first floor( K / 2 ) and last
    % floor( ( K - 1 ) / 2 ) rows
    case  's'
      
      % Number of top rows to chop
      nt = floor (  K  /  2  ) ;
      
      % Number of last rows to chop
      nl = floor (  ( K - 1 )  /  2  ) ;
      
      % Chop rows
      c( [ 1 : nt , end - nl + 1 : end ] , co{ : } ) = [ ] ;
    
    % Causal, remove last K - 1 rows
    case  'c'  ,  c( end - K + 2 : end , co{ : } ) = [ ] ;
      
    % Empty string , no action
    case  ''
      
    % Unrecognised input, fire off a warning but otherwise do nothing
    otherwise
      
      warning (  'MAK:makconv:t_unrecognised'  ,  ...
        'makconv: unrecognised value of t: %s'  ,  t  )
      
  end % kernel type
  
  
end % makconv

