
function  [ k , t ] = makpspkern (  par  )
% 
% k = makpspkern
% k = makpspkern (  par  )
% [ k , t ] = makpspkern (  ...  )
% par = makpspkern (  'default'  )
% 
% MET Analysis Kit. Return a double floating point row vector containing a
% convolution kernel in the shape of a postsynaptic potential. This was
% originally defined by Thompson et al (1996) as a causal kernel for spike
% train smoothing. The kernel is normalised so that it sums to 1. Optional
% output t gives the vector of sample times for each sample in k.
% 
% Without arguments, convolution kernel k is returned using default
% parameters. An optional parameter struct par can be provided to change
% the shape, length, and sampling rate of the kernel in k. Alternatively,
% parameter struct par containing default values can be returned if the
% input argument is a string saying 'default'.
% 
% 
% Parameter struct
% 
% par must be a struct with the following fields:
% 
%   .Tgrowth - A scalar numeric value giving the time constant for the
%     growth phase of the kernel, in seconds. Default 0.001s.
%   
%   .Tdecay - Scalar numeric value with the time constant for the decay
%     phase of the kernel, in seconds. Default 0.020s.
%   
%   .fsample - Sampling frequency of the kernel, in Hertz. Default 1000Hz,
%     giving a kernel with millisecond resolution.
%   
%   .length - Length of the kernel in seconds. The number of samples is
%     .length * .fsample, rounded up. Default 0.2s.
% 
%   All parameters must be real number values greater than 0.
% 
% 
% Example
% 
% % Generate poisson spike train of 0's and 1's at a rate of 30 spk/sec
% spk = rand (  1  ,  1e3  )  <=  30 / 1e3 ;
% 
% % Spike times , in seconds
% tim = ( find (  spk  )  -  1 )  /  1e3 ;
% 
% % PSP kernel
% k = makpspkern ;
% 
% % Convolve spike train
% c = conv (  double( spk )  ,  k  ) ;
% 
% % Throw away convolution tail
% c( 1e3 + 1 : end ) = [] ;
% 
% % Plot convolved spike train
% plot (  ( 1 : 1e3 ) / 1e3  ,  c  ,  'LineWidth'  ,  2  )
% hold  on
% 
% % Plot spike raster
% plot (  [ tim ; tim ]  ,  [ 0 ; max( ylim ) / 3 ]  ,  'k'  )
% 
% 
% Reference
% 
% Thompson, K. G., D. P. Hanes, N. P. Bichot and J. D. Schall (1996).
%   "Perceptual and motor processing stages identified in the activity of
%   macaque frontal eye field neurons during visual search." J Neurophysiol
%   76(6): 4040-4055.
% 
% 
% Written by Jackson Smith - March 2019 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Number of input and output args
  narginchk  (  0  ,  1  )
  nargoutchk (  0  ,  2  )
  
  % No input
  if  nargin  ==  0
    
    % Use default parameters
    par = defpar ;
    
  % Request return of default parameter struct
  elseif  iscellstr (  { par }  )  &&  strcmp (  par  ,  'default'  )
    
    % Return default parameter struct
    k = defpar ;
    return
    
  % Invalid parameter struct provided 
  elseif  validpar ( par )
    
    error (  'MAK:makpspkern:par'  ,  [ 'makpspkern: par must be a ' , ...
      'valid parameter struct, see help makpspkern' ]  )
    
  end % check input
  
  
  %%% Make convolution kernel %%%
  
  % Number of samples
  n = ceil (  par.length  *  par.fsample  ) ;
  
  % Time vector , sample points for total length
  t = 1 ./ par.fsample  *  ( 0 : n ) ;
  
  % Evaluate kernel at time points
  k = ( 1  -  exp( - t ./ par.Tgrowth ) )  .*  exp ( - t  ./  par.Tdecay );
  
  % Normalise so that k sums to 1
  k( : ) = k  ./  sum (  k  ) ;
  
  
end % makpspkern


%%% Sub-routines %%%

% Create default parameter struct
function  par = defpar
  
  % Growth time constant , seconds
  par.Tgrowth = 0.001 ;
  
  % Decay time constant , seconds
  par.Tdecay = 0.020 ;
  
  % Sampling frequency , Hertz
  par.fsample = 1000 ;
  
  % Kernel length , seconds
  par.length = 0.2 ;
  
end % defpar


% Validate parameter struct , returns 0 if parameter struct is valid
function  val = validpar ( par )
  
  % Fail by default
  val = true ;
  
  % Default parameter struct
  def = defpar ;
  
  % Field names
  fdef = fieldnames (  def  ) ;
  
  % Must be a scalar struct with correct field names
  if  ~ isstruct (  par  )  ||  ~ isscalar (  par  )  ,  return  ,  end
  
  % Get field names
  fpar = fieldnames (  par  ) ;
  
  % Check field names
  if  numel (  fdef  )  ~=  numel (  fpar  )  ||  ...
      ~ all (  ismember( fpar , fdef )  )
    
    return
    
  end % field names
  
  % All values must be scalar numeric real and greater than zero
  for  F = fpar'  ,  f = F{ 1 } ;
    
    % Value
    x = par.( f ) ;
    
    % Validate value
    if  ~ isscalar (  x  )  ||  ~ isnumeric (  x  )  ||  ...
        ~ isreal (  x  )  ||  x  <=  0
      
      return
      
    end % validate value
    
  end % value check
  
  % Parameter struct is valid
  val = false ;
  
end % validpar

