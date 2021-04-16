
function  a = makalignspks (  p  ,  t  ,  w  )
% 
% a = makalignspks (  p  ,  t  ,  w  )
% 
% MAK Analysis Kit, pre-processing. Returns the spike waveforms in w
% aligned to their peaks, according to the parameters in struct p. This is
% the first step towards spike sorting.
% 
% Waveforms are converted to single floating point numbers in micro-volts.
% The peak is then found within a time window that follows the theshold
% crossing. Whether the peak is a maximum or minimum is determined based on
% whether the electrode threshold t is greater or less than zero. Once the
% peak sample is found, the centre of mass in the local vicinity of the
% peak is computed in order to estimate the time of the true peak. Once
% peaks and offsets are found, each waveform is resampled using spline
% interpolation ; the result is an aligned waveform.
% 
% 
% Input
% 
% p - Struct parameter with fields:
% 
%   .prethr - Number of samples taken prior to threshold crossing.
%   .peakwin - The number of seconds following the threshold crossing in
%     which to search for a peak.
%   .fs - The sampling rate of the waveforms, in Hertz.
%   .comwin - An index vector for computing centre of weight. Indices are
%     relative to the peak sample.
% 
% t - Electrode threshold level vector, in mico-volts.
% 
% w - S x N matrix of spike waveforms. Each waveform has S samples, indexed
%   along rows. N waveforms are indexed across columns. Expected to be
%   integer numeric values that must be converted to micro-volts.
% 
% 
% Output
% 
% a - ( S - j ) x N matrix of peak-aligned spike waveforms. j is the number
%   of samples in p.peakwin seconds, rounded. Again, rows are indexed by
%   sample and columns are indexed by spike.
% 
% 
% See:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Empty input %%%
  
  % Return an empty single floating point numeric array
  if  isempty ( w )
    a = zeros (  0  ,  'single'  ) ;
    return
  end
  
  
  %%% Find peak samples %%%
  
  % How many samples following the threshold crossing should be checked?
  jit = round ( p.peakwin  *  p.fs ) ;
  
  % Row indices of samples to test for peak
  j = p.prethr : p.prethr + jit ;
  
  % Find maxima following positive threshold crossing
  if  0  <  t
    
    [ ~ , ps ] = max (  w( j , : )  ,  []  ,  1  ) ;
  
  % Find minima following negative threshold crossing
  else
    
    [ ~ , ps ] = min (  w( j , : )  ,  []  ,  1  ) ;
  
  end % Peak sample vector
  
  % Re-align peak row indices to waveforms
  ps = ps  +  p.prethr  -  1 ;
  
  
  %%% Interpolation %%%
  
  % Convert waveforms to single precision floating point numbers, in micro-
  % volts. Pack as a cell vector for ease of data distribution in parfor
  % loop.
  w = single (  w  )  *  p.coef_int2uv ;
  
  % Number of samples in aligned waveform
  na = size (  w  ,  1  )  -  jit ;
  
  % Number of samples to drop from head of raw waveforms
  h = ps  -  p.prethr ;
  
  % Raw waveform coordinate vector
  x = ( 1 : size( w , 1 ) )' ;
  
  % Point to comwin for parfor loop so that whole of p is not copied to
  % each worker
  comwin = p.comwin ;
  
  % Expand into a set of coordinate column vectors
  parfor  i = 1 : numel( h )
    
    % Row indices of waveform samples for cog calculation
    r = ps( i )  +  comwin ;
    
    % Copy single spike waveform to temporary variable
    v = w( : , i ) ;
    
    % Find centre of gravity offset for waveform around its peak sample
    c = comwin  *  v( r )  /  sum (  v( r )  ) ;
    
    % Calculate interpolant coordinates
    xq = ( h( i ) + 1 : h( i ) + na )  +  c ;
    
    % Interpolate
    a( : , i ) = interp1 (  x  ,  v  ,  xq  ,  'spline'  ) ;
    
  end % expand
  
  
end % makalignspks

