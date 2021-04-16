
function  [ Ibc , varargout ] = ...
  maklinfin (  T  ,  N  ,  ds  ,  f  ,  C  ,  fun  )
% 
% Ibc = maklinfin (  T  ,  N  ,  ds  ,  f  ,  C  ,  fun  )
% 
% MET Analysis Kit. Computes bias-corrected linear Fisher information
% according to the method of Kanitscheider et al. ( 2015a ; 2015b ).
% fun is an optional input containing a string that names which sub-
% function to use. Linear Fisher information quantifies how much
% information can be extracted from the population responses by the optimal
% linear decoder. Here, this is done in the context of a stimulus
% discrimination, where stimuli can have values of s = s0 + ds/2 or
% s = s0 - ds/2. For each stimulus value, the responses of N units e.g.
% neurones are recorded over T trials. The average response to each
% stimulus value and the covariance of responses to each value are used to
% compute information, which has the same unit as the stimulus parameter to
% the power of -2. Hence, an orientation discrimination task would return
% linear Fisher information in units of rad ^ -2.
%
% 
% Input arguments
% 
% T , N - Scalar non-zero floating point. The number of trials (T) and the
%   number of units e.g. neurones (N). To be valid, T and N must satisfy
%   this expression: T >= ( N + 5 ) / 2. Otherwise an error is returned.
% 
% ds - Scalar floating point. The difference in stimulus values. Must be a
%   finite value ( see isfinite ). 
% 
% f - Cell array. Each cell contains an N-element vector containing the
%   empirical tuning curve for each stimulus condition. This might be the
%   average firing rate. The exact shape and contents of f depends on fun
%   (see Functions section below). Tuning curves are floating point values.
% 
% C - Cell array. Same size as f. Each cell contains the N x N covariance
%   matrix for each stimulus condition. Covariance matrices are floating
%   point values.
% 
% fun - String. Names what kind of linear Fisher information to compute,
%   and determines the shape of f and C. Valid strings are listed in the
%   Functions section.
% 
% 
% Functions
% 
% These are valid strings for fun.
% 
% 'normal' (default) - The empirical bias-corrected linear Fisher
%   information of the raw data using the optimal linear decoder. This is
%   the default action when fun is omitted. It assumes that population
%   responses were obtained for two separate stimulus values under
%   otherwise identical experimental conditions. This is the Ibc value of
%   equation 2 in Kanitscheider et al. 2015a . f and C are 1 x 2 vectors
%   containing the average firing rate and covariance matrix for each
%   stimulus value.
% 
% 'cross' - The linear Fisher information available when the population
%   responses to two separate stimulus values were obtained in two
%   different experimental conditions, yielding four sets of responses. If
%   a linear decoder is trained on the responses from one experimental
%   condition A then the Fisher information measures the information that
%   can extracted by that decoder when evaluating responses from condtion
%   B. The term Ibc,AB of equation 28 (Kanitscheider et al. 2015a) is thus
%   returned. f and C are 2 x 2 arrays. Row 1 contains data from condition
%   A, and row 2 contains data from condition B. Likewise, column 1
%   contains data for stimulus value 1, and column 2 contains data for
%   stimulus value 2.
% 
% 'shuffle' - The linear Fisher information for one experimental condition
%   and two stimulus values when the responses of each unit are
%   independently shuffled within each stimulus condition. In other words,
%   the amount of information that is encoded by the shuffled responses.
%   Returns the Ibc,shuffle term of equation 6 (Kanitscheider et al.
%   2015a). f and C are the same as for 'normal'.
% 
% 'diag' - The linear Fisher information for one experimental condition and
%   two stimulus values when a linear decoder is optimised for the un-
%   correlated responses but then evaluates the raw, correlated responses.
%   f and is the same as for 'normal'. C has the same format as for
%   'cross', where row 1 contains the shuffled covariance matrices, and row
%   2 has the raw covariance matrices with correlations intact. Returns the
%   Ibc,AB term of equation 29 (Kanitscheider et al. 2015a).
%   
%   Shuffling destroys the correlations between units while maintaining all
%   marginal properties of the data. In the authors' implementation
%   (BCFisherDiag.m), the shuffled covariance matrix is estimated using a
%   single random shuffling of the data. Ideally, the shuffled covariance
%   matrix could be derived from the raw covariance matrix directly, but
%   this would require a new derivation of the bias to be subtracted.
%   Presently, the user must shuffle the data (just once) and recompute the
%   shuffled covariance matrices, to be handed in via argument C.
%   
%   TO DO:
%   
%   In principal, it is preferrable to repeat this process many times with
%   a kind of Monte Carlo method to build up a distribution of values from
%   which the shuffled covariance matrix can be estimated. Otherwise, the
%   estimated shuffled covariance matrix will be noisy. In principal, the
%   off-diagonal elements of the shuffled covariance matrix should approach
%   zero in the limit, while the diagonal elements should match the
%   empirical covariance matrix. A re-derivation of Ibc,diag will be
%   required for this.
% 
% 
% Output arguments
%
% All outputs are scalar floating point.
% 
% Ibc -  The linear Fisher information computed
%   according to the given function in fun.
% 
% [ ... , varIbc , I ] = maklinfin (  ...  ,  'normal'  ) - It is possible
%   to estimate the variance of the Ibc value. This is returned in an
%   optional second output argument, varIbc. If T = ( N + 5 ) / 2 then
%   varIbc will be Inf. See equation 19 (Kanitscheider et al. 2015a). I is
%   the uncorrected linear Fisher information.
% 
% [ ... , varE , biasE ] = maklinfin (  ...  ,  'cross'  ) - Can return
%   varE, the decoder variance on the test data ( wD * SigmaE * wD' ), and
%   also biasE, the decoder bias derivative on the test data 
%   ( wD * muprimeE ).
% 
% [ ... , I ] = maklinfin (  ...  ,  'shuffle'  ) - Also returns the un-
%   corrected linear Fisher information.
% 
%
% References
% 
% Kanitscheider, I., et al. (2015a). "Measuring Fisher Information
%   Accurately in Correlated Neural Populations." PLoS Comput Biol 11(6):
%   e1004218.
% 
% Kanitscheider, I., et al. (2015b). "MatLab tools for estimating linear
%   Fisher information from population data along with synthetic data and
%   recorded spike count responses from neurons in macaque primary visual
%   cortex to grating images with different orientations and white noise".
%   CRCNS.org. http://dx.doi.org/10.6080/K0PK0D3B
%
% 
% Written by Jackson Smith - August 2018 - DPAG , University of Oxford
% 
  

  %%% CONSTANTS %%%
  
  % Function strings
  FUNSTR = {  'normal'  ,  'cross'  ,  'shuffle'  ,  'diag'  } ;
  
  % Max number of output args , per function
  ARGOUT = [        3   ,       3   ,         2   ,      1   ] ;
  
  
  %%% Check input %%%
  
  % Max number of inputs
   narginchk (  5  ,  6  )
  
  % Scalar integer arguments
  for  A = {  { T , 'T' }  ,  { N , 'N' }  }  ,  a = A{ 1 }{ 1 } ;
    
    % Check scalar
    if  ~ isscalar (  a  )  ||  ~ isfloat (  a  )  ||  a  <  0
      
      % Message identifier
      MSGID = sprintf (  'MAK:maklinfin:%s'  ,  A{ 1 }{ 2 }  ) ;
      
      % Launch error
      error (  MSGID  ,  [ 'maklinfin: %s must be scalar floating ' , ...
        'point above zero' ]  ,  A{ 1 }{ 2 }  )
      
    end % check scalar
    
  end % scalar ints
    
  % ds difference in two stimulus values
  if  ~ isscalar (  ds  )  ||  ~ isfloat (  ds  )  ||  ...
      ~ isfinite (  ds  )
    
    error (  'MAK:maklinfin:ds'  ,  ...
      'maklinfin: ds must be scalar, finite, floating point'  )
    
  % Do T and N satisfy equation T >= ( N + 5 ) / 2?
  elseif  T  <  ( N + 5 )  /  2
    
    error (  'MAK:maklinfin:TN'  ,  ...
      'maklinfin: T and N do not satisfy expression T >= ( N + 5 ) / 2'  )
    
  end % scalar arg checks
  
  % Cell array input args
  for  A = {  { f , @any , 'N-element vector' , 'f' }  ,  ...
              { C , @all , 'N x N matrix'     , 'C' }  }
    
    % Argument value
    a = A{ 1 }{ 1 } ;
    
    % Message identifier
    MSGID = sprintf (  'MAK:maklinfin:%s'  ,  A{ 1 }{ 4 }  ) ;
    
    % No more than 2 dims allowed
    if  ~ ismatrix (  a  )
      
      error (  MSGID  ,  [ 'maklinfin: %s must have no more than ' , ...
          '2 dimensions' ]  ,  A{ 1 }{ 4 }  )
      
    end % 2 dims
    
    % Each cell
    for  i = 1 : numel (  a  )
      
      % Size of contents
      s = size (  a{ i }  ) ;
      
      % Check
      if  ~ isfloat (  a{ i }  )  ||  2  <  numel (  s  )  ||  ...
          ~ A{ 1 }{ 2 }(  s  ==  N  )

        % Launch error
        error (  MSGID  ,  [ 'maklinfin: %s{ %d } must be a ' , ...
          'floating point %s' ]  ,  A{ 1 }{ 4 }  ,  i  ,  A{ 1 }{ 3 }  )
        
      end % contents check
      
    end % cells
    
  end % cell arrays
  
  % Function string given?
  if  5  <  nargin
    
    % Is it a string?
    if  ~ iscellstr (  { fun }  )
      
      error (  'MAK:maklinfin:fun_str'  ,  ...
        'maklinfin: fun must be a string'  )
      
    end % fun is string
    
    % Find matching function string
    i = strcmp(  fun  ,  FUNSTR  ) ;
    
    % Valid function string?
    if  ~ any (  i  )
      
      error (  'MAK:maklinfin:fun_val'  ,  ...
        'maklinfin: fun has unrecognised function string: %s'  ,  fun  )
      
    end % valid fun str
    
  % No!
  else
    
    % Use default
    fun = 'normal' ;
    
    % Find matching function string
    i = strcmp(  fun  ,  FUNSTR  ) ;
    
  end % function string
  
  % Check max number of output arguments
  nargoutchk (  0  ,  ARGOUT( i )  )
  
  % Get required size of f and C according to function used
  switch  fun
    
    % We require 2 x 2 cell arrays for f and C input to 'cross'
    case  'cross'
      
      ROW.f = 2 ;  COL.f = 2 ;
      ROW.C = 2 ;  COL.C = 2 ;
      
    % We require 2 x 2 cell arrays for C but not f for input to 'diag'
    case  'diag'
      
      ROW.f = 1 ;  COL.f = 2 ;
      ROW.C = 2 ;  COL.C = 2 ;
      
    % All other cases require 1 x 2 cell arrays
    otherwise
      
      ROW.f = 1 ;  COL.f = 2 ;
      ROW.C = 1 ;  COL.C = 2 ;
    
  end % cell array input size
  
  % Check sizes
  chksize (  fun  ,  f  ,  'f'  ,  ROW.f  ,  COL.f  )
  chksize (  fun  ,  C  ,  'C'  ,  ROW.C  ,  COL.C  )
  
  
  %%% Compute linear Fisher information %%%
  
  % Sub-functions
  switch  fun
    
    
    % Empirical Fisher info from optimal linear decoder for raw responses
    case  'normal'
      
      % Average covariance matrix between two stimulus conditions
      S = (  C{ 1 }  +  C{ 2 }  )  /  2 ;
      
      % Difference in empirical tuning curves between stim conditions
      df = (  f{ 1 }( : )  -  f{ 2 }( : )  )'  /  ds ;
      
      % Naive linear Fisher information , no bias correction
      I = df  *  pinv (  S  )  *  df' ;
      
      % Bias-corrected linear Fisher information
      Ibc = I  *  (  ( 2 * T - N - 3  )  /  ( 2 * T - 2 )  )  -  ...
        ( 2 * N )  /  ( T * ds ^ 2 ) ;
      
      % Variance estimator , if requested
      if  1  <  nargout
        
        varargout{ 1 } = ( 2 * Ibc ^ 2 )  /  ( 2 * T - N - 5 )  *  ...
          ( 1  +  4 * ( 2 * T - 3 ) / ( T * Ibc * ds ^ 2 ) + ...
            4 * N * ( 2 * T - 3 ) / ( T * Ibc * ds ^ 2 ) ^ 2 ) ;
        
      end % VarIbc
      
      % Naive linear Fisher info , if requested
      if  2  <  nargout  ,  varargout{ 2 } = I ;  end
      
      
    % Fisher info of linear decoder trained in condition A and evaluating
    % responses from condition B
    case  'cross'
      
      % Difference of tuning curves for decoder training condition
      emupA = (  f{ 1 , 1 }( : )  -  f{ 1 , 2 }( : )  )'  /  ds ;
      
      % Difference of tuning curves for decoder testing condition
      emupB = (  f{ 2 , 1 }( : )  -  f{ 2 , 2 }( : )  )'  /  ds ;
      
      % Decoder training condition covariance matrix
      eCA = (  C{ 1 , 1 }  +  C{ 1 , 2 }  )  /  2 ;
      
      % Decoder testing condition covariance matrix
      eCB = (  C{ 2 , 1 }  +  C{ 2 , 2 }  )  /  2 ;
      
      % Bias-corrected inverse covariance matrix for decoder training
      % condition
      invCA = pinv (  eCA  )  *  ( 2 * T - N - 3 )  /  ( 2 * T - 2 ) ;
      
      % Bias-corrected information for optimal decoder of responses in
      % training condition
      FIA = ( emupA * invCA * emupA' )  -  2 * N / ( T * ds ^ 2 ) ;
      
      % Prepare terms for computing lambda and rho
      d = ( 2 * T - N - 2 )  *  ( 2 * T - N - 5 ) ;
      
      tr = trace (  eCB  *  invCA  ) ;
      
      % Lambda and rho , see equations 27 of Kanitscheider et al. (2015a)
      lambda = 2  /  ( T * ds ^ 2 )  *  ( tr * ( 1 + ( ...
        ( 2 * T - N - 1 ) + N * ( 2 * T - N - 3 ) ) / d ) ) ;
      
      rho = FIA  *  tr  *  ( 2 * T - N - 3 )  /  d ;
      
      % Un-corrected linear Fisher info numerator
      uFInum = emupA * invCA * emupB' ;
      
      % Term used to compute both the corrected Fisher info and the
      % variance
      VARnum = (  ( emupA * invCA * eCB * invCA * emupA' )  -  ...
        lambda  -  rho  )  /  (  ( 1 + ( 2 * T - N - 1 ) / d )  ) ;
      
      % Bias-corrected estimate of the crossed information in the testing
      % responses using the now suboptimal decoder. Output FIE in
      % BCFisherCross.m
      Ibc =   uFInum ^ 2 / VARnum ;
      
      % Variance estimator , if requested
      if  1  <  nargout  ,  varargout{ 1 } = VARnum  /  FIA ^ 2 ;  end
      
      % Bias estimator , if requested
      if  2  <  nargout  ,  varargout{ 2 } = uFInum  /  FIA     ;  end
      
      
    % Fisher info of linear decoder optimised for shuffled responses and
    % evaluated on shuffled responses
    case  'shuffle'
      
      % Average variance between two stimulus conditions
      V = (  diag(  C{ 1 }  )  +  diag(  C{ 2 }  )  )  /  2 ;
      
      % Difference in empirical tuning curves between stim conditions. Take
      % squared differences divided by squared difference in stim values.
      df = (  f{ 1 }( : )  -  f{ 2 }( : )  )  .^  2  /  ds ^ 2 ;
      
      % Un-corrected Fisher info for each unit
      Ieach = df  ./  V ;
      
      % Eliminate units with non-finite info values e.g. when no spikes
      % were fired and variance was zero
      Ieach(  ~ isfinite(  Ieach  )  ) = [] ;
      
      % How many neurones remain?
      N = numel (  Ieach  ) ;
      
      % Un-corrected population shuffled linear Fisher info
      I = sum (  Ieach  ) ;
      
      % Bias-corrected population shuffled linear Fisher info
      Ibc = I  *  ( T - 2 )  /  ( T - 1 )  -  2  *  N  /  ( T * ds ^ 2 ) ;
      
      % Naive linear Fisher info , if requested
      if  1  <  nargout  ,  varargout{ 1 } = I ;  end
      
      
    % Fisher info of linear decoder optimised for shuffled responses but
    % evaluated on raw responses with correlations intact
    case  'diag'
      
      % Average covariance matrix between two stimulus conditions ,
      % shuffled data dampens correlations . This is the training condition
      % for the crossed information. Hence, condition A.
      eCA = (  C{ 1 , 1 }  +  C{ 1 , 2 }  )  /  2 ;
      
      % Average covariance matrix between two stimulus conditions , raw
      % data and correlations are intact. This is the testing condition for
      % the crossed information. Hence, condition B.
      eCB = (  C{ 2 , 1 }  +  C{ 2 , 2 }  )  /  2 ;
      
      % Difference in empirical tuning curves between stim conditions
      df = (  f{ 1 }( : )  -  f{ 2 }( : )  )'  /  ds ;
      
      % [ If Ibc,shuffle is rederived properly then remove the eCA above
      % and uncomment this section ]
      %
%       % Shuffled covariance matrix , un-correlated data. This is the
%       % training condition for crossed information. Hence, condition A.
%       eCA = zeros ( size(  eCB  ) ) ;
%       
%       % Diagonal elements
%       i = find ( eye(  N  ) ) ;
%       
%       % Transfer variances
%       eCA( i ) = eCB( i ) ;
      
      % Naive, un-corrected inverse covariance of shuffled data
      einvCA = pinv (  eCA  ) ;
      
      % Bias-corrected inverse covariance of shuffled data
      invCD = einvCA  *  ( 2 * T - N - 3 )  /  ( 2 * T - 2 ) ;
      
      % Shuffled information term , used for next and last steps
      FIAterm = df  *  invCD  *  df' ;

      % Bias-corrected estimate of information for shuffled data (optimal
      % decoder)
      FIA = FIAterm  -  2 * N / ( T * ds ^ 2 ) ;

      % Prepare some simplifying terms that are used in several steps
      c = 2  *  ( T - 1 )  /  ( 2 * T - N - 3 ) ;
      d = ( 2 * T - N - 2 )  *  ( 2 * T - N - 5 ) ;
      tr = trace (  eCB  *  invCD  ) ;
      
      % Prepare simplifying terms used in the next step
      a = ( ( 2 * c ^ 2 ) / ( T * ds ^ 2 ) )  *  ( tr * ( 1 + ( ...
        ( 2 * T - N - 1 ) + N * ( 2 * T - N - 3 ) ) / d ) ) ;
      b = FIA  *  tr  *  ( c ^ 2 )  *  ( 2 * T - N - 3 )  /  d ;
      
      % Wibbly-wobbly term required for the last step
      VARnum = ( ( df * einvCA * eCB * einvCA * df' ) - a - b )  /  ...
        ( c ^ 2 * ( 1 + ( 2 * T - N - 1 ) / d ) ) ;
      
      % Bias-corrected linear Fisher information of correlated responses
      % evaluated by suboptimal decoder trained on un-correlated responses
      Ibc = (  FIAterm  -  2 * N / ( T * ds ^ 2 )  ) ^ 2  /  VARnum ;
      
      
    % Haven't implemented this function yet, though it is listed in set of
    % function names
    otherwise
      
      error (  'MAK:maklinfin:fun_imp'  ,  ...
        'maklinfin: function not implemented fun: %s'  ,  fun  )
      
      
  end % sub-functions
  
  
end % maklinfin


%%% Sub-routines %%%

function  chksize (  fun  ,  a  ,  s  ,  ROW  ,  COL  )
  
  % Get number of rows and columns
  [ r , c ] = size (  a  ) ;

  % Check size of cell array
  if  r  ~=  ROW  ||  c  ~=  COL

    error (  'MAK:maklinfin:cell_size' ,  [ 'maklinfin: function ' , ...
      '''%s'' requires %s to be %d x %d' ] ,  fun ,  s ,  ROW ,  COL  )

  end % cell array size
  
end % chksize

