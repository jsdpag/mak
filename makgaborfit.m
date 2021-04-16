
function  [ C , r2 , lb , ub , c0 ] = makgaborfit(  x  ,  Y  ,  varargin  )
% 
% C = makgaborfit (  x  ,  Y  )
% C = makgaborfit (  x  ,  Y  ,  Ya  )
% [ C , r2 , lb , ub , C0 ] = makgaborfit (  ...  )
% pars = makgaborfit
% ... = makgaborfit( x , Y , ... , pars )
% 
% MET Analysis Kit. Uses non-linear, least-squares optimisation to find the
% 1-dimensional Gabor function that best predicts values in matrix Y from
% values in vector x. The size of dimension 1 in Y must equal to the length
% of x (column major number of elements). Each column of Y is considered to
% be a different data set, so a Gabor is fit to each column of data.
% Coefficients of best-fit Gabors are returned in matrix C, with a column
% for each column of Y.
% 
% Rows of C have coefficients in the order:
%   
%   1) y0 - Baseline value
%   2) A  - Gaussian envelope amplitude
%   3) x0 - Gaussian envelope horizontal offset
%   4) s  - Width of the Gaussian envelope (s for sigma)
%   5) f  - Frequency of cosine
%   6) p  - Phase of cosine , in radians
% 
% If there is an auxiliary, corresponding set of data then it can be passed
% in as the optional third argument, Ya. Ya must have the same size as Y,
% and the data Y( : , i ) and Ya( : , i ) are used in the same optimisation
% procedure to constrain the parameter search. In this case, an additional
% two parameters are appended in two final rows of C:
% 
%   7) Aa - Auxiliary Gaussian envelope amplitude
%   8) pa - Auxiliary phase of cosine, in radians
% 
% Thus, when Ya is given then two Gabors are fit using the same y0, x0, s,
% and f values but different amplitude and phase values i.e. the amplitude
% and phase are varied separately for each data set.
%
% As a bonus, you can also get the coefficient of determination in r2, a
% row vector with a measurement for each fitted Gabor. This is defined as:
% 
%   r2 = 1  -  SSres / SStot
%   
%   where SSres is the residual sum of squares and SStot is the total sum
%   of squares
% 
% If Ya is provided, then r2 becomes a 2 x N array, where row 1 has
% coefficients for Gabors fit to Y and row 2 has coefficients for Ya.
% 
% makgaborfit can return the lower and upper bounds of the parameter search
% in lb and ub, and the starting parameters in C0. These are all the same
% size as C, with the same organisation.
% 
% This procedure is optimised for fitting disparity tuning curves following
% the Methods of:
%
%   Prince, S. J. D., A. D. Pointon, B. G. Cumming and A. J. Parker (2002).
%     Quantitative Analysis of the Responses of V1 Neurons to Horizontal
%     Disparity in Dynamic Random-Dot Stereograms. Journal of
%     Neurophysiology 87(1): 191-208.
%   
%   Tanabe S, Umeda K, Fujita I. Rejection of false matches for binocular
%     correspondence in macaque visual cortical area V4. J Neurosci. 2004
%     Sep 15;24(37):8170-80.
%   
%   Cumming, B. G. and A. J. Parker (1997). Responses of primary visual
%     cortical neurons to binocular disparity without depth perception.
%     Nature 389(6648): 280-283.
% 
% Hence, x is generally assumed to be a vector of tested binocular
% disparities spanning -1 to +1 at regular steps. However, this will still
% work without regular steps and a narrower range e.g. x = [ -0.5, -0.2,
% -0.1, -0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 ].
% 
% makgaborfit will return a parameter struct in pars that can be altered
% and returned as the final input argument in order to customise the
% bounds of the parameter search. Supported fields of pars includes:
% 
%   pars.xlimits_mult - 2-element vector giving multiplication factors on
%     the difference dx = max( x ) - min( x ). The full range of horizontal
%     offsets will be min( x ) - dx * pars.xlimits_mult( 1 ) up to max( x )
%     + dx * pars.xlimits_mult( 2 ). Valid values are 0 and greater.
%     Default [ 0 , 0 ].
%   
%   pars.width_mult - 2-element vector for lower and upper multiplication
%     factors on range of Gaussian widths. These are based on the
%     difference in upper and lower horizontal offset limits dxlim =
%     x_upper - x_lower. dxlim * pars.width_mult( 1 ) is the lower limit on
%     width and dxlim * pars.width_mult( 2 ) is the upper limit on width.
%     Values must be 0 or greater. Default [ 0.1 , 1 ].
%   
%   pars.amplitude_mult - Amplitude multiplication factor. The upper bound
%     on amplitude is this times the difference between the maximum and
%     minimum observed value. Must be scalar value of 0 or greater. Default
%     2.
%   
%   pars.percent_max_freq - 2-element vector. These are percentage values.
%     The dominant 
%     frequency of the data in Y is found and used as the starting value
%     for coefficient f, the cosine frequency. The lower and upper bounds
%     of frequency in the paremeter will be f - f *
%     pars.percent_max_freq( 1 ) and f + f * pars.percent_max_freq( 2 ).
%     Valid range of values is 0 and greater. Default [ 10% , 10% ].
% 
% See also: makgabor
% 
% Written by Jackson Smith
% 
  
  
  %%% Check input %%%
  
  % Number of allowable input/output args
   narginchk( 0 , 4 )
  nargoutchk( 0 , 5 )
  
  % No input args, return parameter struct
  if  ~ nargin  ,  C = defpar ;  return  ,  end
  
  % Number of points in x
  nx = numel (  x  ) ;
  
  if  ~ isvector (  x  )
    
    error (  'MAK:makgaborfit:x' ,  'makgaborfit: x must be a vector'  )
    
  elseif  ~ isnumeric( x )
    
    error (  'MAK:makgaborfit:num_x' ,  'makgaborfit: x must be numeric'  )
    
  elseif  ~ ismatrix (  Y  )
    
    error (  'MAK:makgaborfit:Y' ,  'makgaborfit: Y must be a 2D matrix'  )
    
  elseif  nx  ~=  size (  Y  ,  1  )
    
    error (  'MAK:makgaborfit:size' ,  ...
      'makgaborfit: length of x must equal size of dim 1 of Y'  )
    
  elseif  ~ isnumeric( Y )
    
    error (  'MAK:makgaborfit:num_Y' ,  'makgaborfit: Y must be numeric'  )
    
  end % check input
  
  % Guarantee that x is a column vector
  if  ~ iscolumn (  x  )  ,  x = x( : ) ;  end
  
  % Argument Ya is provided. 3rd argument must be Ya because the only other
  % valid possibility is the parameter struct.
  if  3 <= nargin  &&  ~ isstruct( varargin{ 1 } )
    
    % Raise auxiliary data flag
    auxflg = true ;
    
    % Point to Ya with meaningful name
    Ya = varargin{ 1 } ;
    
    % Check Ya
    if  ~ ismatrix (  Ya  )
      
      error (  'MAK:makgaborfit:Ya' ,  ...
        'makgaborfit: Ya must be a 2D matrix'  )
      
    elseif  ~ isnumeric( Y )
    
      error (  'MAK:makgaborfit:num_Ya' ,  ...
        'makgaborfit: Ya must be numeric'  )
      
    elseif              nx  ~=  size (  Ya  ,  1  )  ||  ...
        size (  Y  ,  2  )  ~=  size (  Ya  ,  2  )
      
      error (  'MAK:makgaborfit:Ya_size' ,  ...
        'makgaborfit: Ya must be same size as Y'  )
      
    end % chk Ya
    
  % Ya not given
  else
    
    % Lower auxiliary data flag
    auxflg = false ;
    
  end % arg Ya
  
  % Get parameter struct
  if  ~ isempty( varargin )  &&  isstruct( varargin{ end } )
    
    pars = varargin{ end } ;
    
    if  ~ isscalar( pars )
      
      error (  'MAK:makgaborfit:pars_size' ,  ...
        'makgaborfit: pars must be a scalar struct'  )
      
    elseif  ~ all(  isfield( pars , { 'xlimits_mult' , ...
        'width_mult' , 'amplitude_mult' , 'percent_max_freq' } )  )
      
      error (  'MAK:makgaborfit:pars_missing' ,  ...
        'makgaborfit: struct pars is missing fields'  )
      
    end
    
  % Default parameters
  else
    
    pars = defpar ;
    
  end % params
  
  % Check horizontal offset multipliers
  if  ~ validpar( pars.xlimits_mult , 2 )  ||  any( pars.xlimits_mult < 0 )
    
    error (  'MAK:makgaborfit:xlimits_mult' ,  [ 'makgaborfit: ' , ...
      'xlimits_mult are not a valid horizontal offset factors' ]  )
      
  end % xlimits_mult
  
  % Check horizontal offset multipliers
  if  ~ validpar( pars.width_mult , 2 )  ||  any( pars.width_mult < 0 )
    
    error (  'MAK:makgaborfit:width_mult' ,  [ 'makgaborfit: ' , ...
      'width_mult are not a valid width factors' ]  )
      
  end % xlimits_mult
  
  % Check amplitude_mult
  if  ~ validpar( pars.amplitude_mult , 1 )  ||  pars.amplitude_mult < 0
    
    error (  'MAK:makgaborfit:amplitude_mult' ,  ...
        'makgaborfit: amplitude_mult is not a valid amplification'  )
    
  end % amplitude_mult
  
  % Check percent_max_freq is a percentage value. If so then convert to
  % fraction.
  if  validpar( pars.percent_max_freq , 2 )  &&  ...
      all( 0 <= pars.percent_max_freq )
    
    pars.percent_max_freq = pars.percent_max_freq  ./  100 ;
    
  else
    
    error (  'MAK:makgaborfit:percent_max_freq' ,  ...
       'makgaborfit: percent_max_freq are not a valid percentage scores'  )
    
  end % percent_max_freq
  
  
  %%% Constants %%%
  
  % lsqcurvefit options object , needed to disable messages
  opt = optimoptions (  'lsqcurvefit'  ,  'Display'  ,  'none'  ) ;
  
  % Number of coefficients
  Nc = 6 ;
  
    % Auxiliary data Ya given , allocate additional amp and phase terms
    if  auxflg  ,  Nc = Nc  +  2 ;  end
  
  % Number of data sets
  Nd = size (  Y  ,  2  ) ;
  
  % Number of interpolant points
  Ni = 2  ^  10 ;
  
  % Coefficient of determination was requested
  if  nargout  >  1
    
    % Raise flag
    r2flg = true ;
    
  else
    
    % Lower flag
    r2flg = false ;
    
  end % r2 wanted
  
  
  %%% Allocate output %%%
  
  C = zeros (  Nc  ,  Nd  ) ;
  
  % No data to fit , end here
  if  ~ Nd  ,  return  ,  end
  
  
  %%% Coefficient limits and starting values %%%
  
  % Allocate lower and upper bound matrices , and starting coefficients
  lb = zeros (  Nc  ,  Nd  ) ;
  ub = zeros (  Nc  ,  Nd  ) ;
  c0 = zeros (  Nc  ,  Nd  ) ;
  
  % 1) Baseline offset constrained between zero and the maximum observed
  %    response
  ub( 1 , : ) = max (  Y  ,  [ ]  ,  1  ) ;
  
    % Use average across all values of x as starting point
    c0( 1 , : ) = mean (  Y  ,  1  ) ;
  
  % 2) Amplitude constrained between zero and twice the difference between
  %    the observed maximum and minimum
  
    % Use that difference as the start point
    c0( 2 , : ) = ub( 1 , : )  -  min (  Y  ,  [ ]  ,  1  ) ;
    
    % Then compute the upper bound
    ub( 2 , : ) = pars.amplitude_mult * (  c0( 2 , : )  ) ;
  
  % 3) Minimum and maximum horizontal offsets are constrained to the range
  %    of values of x. Starting values is median of x.
  xmin = min (  x  ) ;
  xmax = max (  x  ) ;
    dx = xmax  -  xmin ;
  
    % Bounds
    lb( 3 , : ) = xmin  -  dx * pars.xlimits_mult( 1 ) ;
    ub( 3 , : ) = xmax  +  dx * pars.xlimits_mult( 2 ) ;
    
      % Difference in bounds
      dx = ub( 3 , : )  -  lb( 3 , : ) ;
    
    % Starting value
    c0( 3 , : ) = median (  x  ) ;
  
  % 4) Width of Gaussian envelope constrained between 0.1 and total range
  %    of values in x
  lb( 4 , : ) = pars.width_mult( 1 )  *  dx ;
  ub( 4 , : ) = pars.width_mult( 2 )  *  dx ;
  
    % Starting values are half the total range
    c0( 4 , : ) = ub( 4 , : )  /  2 ;
  
  % 5) Estimate frequency in x domain from empirical data and constrain
  %    search to +/- 10% of that.
  
    % Increment per dample
    dx = ( xmax - xmin )  ./  ( Ni - 1 ) ;
    
    % Sampling frequency
    fs = 1  ./  dx ;
  
    % Interpolant values of x
    dint = ( xmin : dx : xmax )' ;
    
    % Take z-score transformation of data
    V = zscore (  Y  ,  0  ,  1  );
    
    % Linear interpolation of data
    VQ = interp1 (  x  ,  V  ,  dint  ) ;
    
    % Fourier transform
    y = fft (  VQ  ) ;
    
    % Spectral power
    p = ( abs( y ) .^ 2 )  /  Ni ;
    
    % Frequency components
    f = ( 0 : Ni - 1 )'  *  ( fs / Ni ) ;
    
    % Maximum frequency component index from 0 to Nyquist of original
    % maximum sampling rate
    [ ~ , j ] = max (  p( 0 < f & f <= 50 , : )  ) ;
    
    % Peak frequency for each data set, and a percentage of that
    df = f( j + 1 ) ;
    
    % Constrain to +/- a percentage of peak frequency
    lb( 5 , : ) = df  -  pars.percent_max_freq( 1 ) * df ;
    ub( 5 , : ) = df  +  pars.percent_max_freq( 2 ) * df ;
    
    % Use disparity frequency as starting point
    c0( 5 , : ) = df ;
    
  % 6) Cosine phase constrained to +/- 3pi. Starting values all zero.
  lb( 6 , : ) = - 3 * pi ;
  ub( 6 , : ) = + 3 * pi ;
  
  % Auxiliary data given
  if  auxflg
    
    % Copy amplitude and phase bounds and starting values
    lb( 7 : 8 , : ) = lb( [ 2 , 6 ] , : ) ;
    ub( 7 : 8 , : ) = ub( [ 2 , 6 ] , : ) ;
    c0( 7 : 8 , : ) = c0( [ 2 , 6 ] , : ) ;
    
  end % Ya
  
  
  %%% Estimate coefficients %%%
  
  % Auxiliary data Ya given
  if  auxflg
    
    % Row indices for auxiliary data, placing auxiliary amp and phase in
    % correct order amongst shared parameters
    riaux = [ 1 , 7 , 3 : 5 , 8 ] ;
    
    % Function evaluates Gabor for each set of parameters
    fg = @( c , x ) [  makgabor( c( 1 : 6 ) , x )  ;
                       makgabor( c( riaux ) , x )  ] ;
                     
    % Append Y and Ya for use by curve fitting function
    Y = [  Y  ;  Ya  ] ;
    
  % No auxiliary data
  else
    
    % No changes to Y, and pass in basic gabor function
    fg = @makgabor ;
    
  end % Ya
  
  % r2 wanted, allocate space for evaluating fitted gabors
  if  r2flg  ,  Yhat = zeros ( nx , Nd , 1 + auxflg ) ;  end
  
  % Data sets
  for  i = 1 : Nd
    
    % Non-linear, least-squares fitting
    C( : , i ) = lsqcurvefit (  fg ,  c0( : , i ) ,  x ,  Y( : , i ) ,  ...
      lb( : , i ) ,  ub( : , i ) ,  opt  ) ;
    
    % r2 wanted
    if  r2flg
      
      % Evaluate Gabor fitted to input data Y
      Yhat( : , i , 1 ) = makgabor (  C( 1 : 6 , i )  ,  x  ) ;
      
      % Auxiliary data given
      if  auxflg
        
        % Evaluate Gabor fitted to auxiliary data Ya
        Yhat( : , i , 2 ) = makgabor (  C( riaux , i )  ,  x  ) ;
        
      end % aux data
    
    end % r2 wanted
   
  end % data sets
  
  % r2 not wanted, quit now
  if  ~ r2flg
    
    return
    
  % But if Ya given then re-shape Y so that the original Y and Ya arguments
  % are appended across dim 3, rather than dim 1 as they are at this point
  elseif  auxflg
  
    Y = cat (  3  ,  Y( 1 : nx , : )  ,  Ya  ) ;
    
  end % proceed to calculating r2
  
  % Residual sum of squares
  SSres = sum (  ( Y - Yhat ) .^ 2  ,  1  ) ;
  
  % Total sum of squares, remember that starting baseline is the mean
  SStot = sum (  bsxfun( @minus , Y , c0( 1 , : ) ) .^ 2  ,  1  ) ;
  
  % Coefficient of determination
  r2 = 1  -  SSres ./ SStot ;
  
  % Auxiliary data given
  if  auxflg
    
    % Reshape r2 into 2 x N matrix, rather than the 1 x N x 2 array that it
    % is
    r2 = permute (  r2  ,  [ 3 , 2 , 1 ]  ) ;
    
  end % Ya
  
end % makgaborfit


%%% Sub-routines %%%

% Define default bounds for parameter search
function  pars = defpar
  
  % Horizontal offset factors
  pars.xlimits_mult = [ 0 , 0 ] ;
  
  % Width multiplication factor
  pars.width_mult = [ 0.1 , 1 ] ;
  
  % Amplitude multiplication factor. Default 2.
  pars.amplitude_mult = 2 ;
  
  % ±percentage of starting frequency for cosine. Default 10%.
  pars.percent_max_freq = [ 10 , 10 ] ;
  
end % defpar

% Check that all elements of x are valid parameters. Scalar, finite value,
% and numeric. chksca if non-zero then checks that x has chknum elements
function  p = validpar( x , chknum )
  
  if  nargin > 1  &&  chknum
    numflg = numel( x ) == chknum ;
  else
    numflg = true ;
  end
  
  p = numflg  &&  all(  isfinite( x )  &  isnumeric( x )  ) ;
  
end % validpar

