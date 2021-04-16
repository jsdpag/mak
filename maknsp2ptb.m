
function  [ n2pcoef , n2pvolts , n2pclim , n2perr ] = ...
            maknsp2ptb (  MSID  ,  epoch  ,  par  ,  nev  ,  nsx  ,  ptb  )
% 
% [ n2pcoef , n2pvolts , n2pclim , n2perr ] =
%           maknsp2ptb (  MSID  ,  epoch  ,  par  ,  nev  ,  nsx  ,  ptb  )
% 
% MET Analysis Kit, pre-processing. Performs robust linear regression
% on Cerebus NSP frame onset times measured with a photodiode to predict
% reported PTB frame onset times. Coefficients are returned in n2pcoef as a
% 2-element row vector ; n2pcoef( 1 ) is the y-intercept and n2pcoef( 2 )
% is the slope. Thus, n2pcoef can be applied to any Cerebus NSP time stamps
% (once converted from samples to seconds) to find out when those events
% happened according to the MET/PTB clock.
% 
% The approach is to use the labelled TTL pulses to identify the end of the
% trial. Following this, photodiode minima and maxima are measured. The
% resulting distributions are used to identify frames that were recorded
% before the start of the trial. The initial synchronising flip is hence
% the first one to exceed the classification limit of pre-trial frames.
% Frames that followed the sync flip are then timed up to the final
% stimulus frame of the trial ; again, these are identified in comparison
% with peak measurements following the trial.
% 
% Once measured frame onset times are obtained, they are regressed against
% the PTB frame onset times i.e. stimulus onset times using robustfit.
% 
% If there were skipped frames then consecutive frames without a skip are
% searched for in two blocks, one from the start of the trial going
% forwards in time, and another from the end of the trial going backwards.
% The number of non-skip frames found in this way must exceed some minimum
% count. If so, then the regression is performed as normal with the
% available non-skip frames.
% 
% However, if there are too many frame skips near the start and end of the
% trial, then the analysis is aborted. mcalibrate time stamps may
% alternatively be used.
% 
% n2pvolts returns the set of in-trial peaks while n2pclim contains the
% classification limits, in raw integer units. Convert both to floating
% point numbers and multiply by p.spksort.coef_int2uv where p is a
% parameter struct from p = makprep ;
% 
% 
% Input arguments
% 
% MSID - A struct that maps MET signal names to their signal identifiers.
% 
% epoch - A scalar logical. Non-zero when the analysis epoch was captured
%   by the the NSP recording.
% 
% par - Parameter struct with at least the following fields:
% 
%   .eid - Electrode ID from .ns* data that carried photodiode
%     measurements.
%   .sync - String saying if the synchronising flip is a maximum ('max') or
%     minimum ('min'). 
%   .flips - String that says whether in-trial flips are identified with
%     the maximum ('max') or minimum ('min') of the photodiode waveform.
%   .leadpeak - If true then time each flip on the peak that immediately
%     preceeds peaks found by .flips.
%   .tailburnmax - Number of peaks to skip following end of trial mstop
%     signal before sampling between-trial mid-grey maxima.
%   .tailburnmin - Same as .tailburnmax, but for minima.
%   .useheader - Instead of estimating the inter-trial peak values from the
%     end of the trial, look at peaks before the trial starts. The
%     .tailburn* parameters will still apply, but in the opposite
%     direction. Set true to use header, defaults false to use tail.
%   .relaxedmatch - If we think that the last frame or two will be lost to
%     the photodiode recording for a known systematic reason, then make
%     this true and maknsp2ptb will match the first N peaks between the
%     recording and the PTB frame times. Default false.
%   .minprom - Minimum allowable prominence of peaks found by findpeaks, in
%     raw units.
%   .minwid - Minimum allowable width of peaks found by findpeaks, in
%     seconds.
%   .stdevs - The number of standard deviations above and below the median
%     post-trial peaks used to define classification boundaries.
%   .skips - If non-zero then trials with frame skips are used.
%   .skipn - The minimum number of non-skip frames that are required to
%     estimate the regression.
%   .usecal - There are too few non-skip frames for the regression, but
%     mcalibrate time stamps may be used.
% 
% nev - Struct returned by Blackrock Microsystem's NPMK function openNEV.
%   Contains trial event data.
% 
% nsx - Struct returned by BM's NPMK function openNSx. Contains photodiode
%   measurements.
% 
% ptb - Struct containing PTB frame times. Variables loaded from the
%   ptbframes_*.mat file.
% 
% 
% n2perr is zero unless an error was detected. On error, n2perr is returned
% with the following values:
% 
%   1 - Missing data. This occurs if either the mstart signal was not
%       recorded by the Cerebus NSP, or if the first part of the photodiode
%       measurements is missing.
%   
%   2 - There are too few non-skip frames to compute the regression, and
%       mcalibrate times are not being used. Or there are too few tail-end
%       off-trial peaks.
%   
%   3 - Unable to perform the regression but mcalibrate times may be used.
% 
% 
% NOTE: Currently assumes that a maxima always preceeds a minima for each
%   frame when par.leadpeak is non-zero. A better way might be to sort all
%   peak times first.
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Constants %%%
  
  % Error saying that there is missing data
  EMISSING = 1 ;
  
  % Error saying that regression was not possible. Returns 2 if mcalibrate
  % is not in use, or 3 if it is.
  ENOREGRESS = 2  +  par.usecal ;
  
  % This flag says that the analysis epoch was recorded and that skipped
  % frames are allowed
  EPOCHF = epoch  &&  par.skips ;
  
  
  %%% Default output arguments %%%
  
  n2pcoef = [ ] ;
  n2pvolts = [ ] ;
  n2pclim = [ ] ;
  n2perr  =  0  ;
  
  
  %%% mstart and mstop signals on NSP clock %%%
  
  % Locate the mstart signal's time
  i = nev.Data.SerialDigitalIO.UnparsedData  ==  MSID.mstart ;
  mstart = nev.Data.SerialDigitalIO.TimeStampSec(  i  ) ;
  
  % Locate the mstop signal's time
  i = nev.Data.SerialDigitalIO.UnparsedData  ==  MSID.mstop ;
  mstop  = nev.Data.SerialDigitalIO.TimeStampSec(  i  ) ;
  
  % It is possible that the Cerebus somehow stores the mstop from the last
  % trial into the recording for this one , check for this in both mst*
  if  ~ isscalar (  mstart  )  ,  mstart = max (  mstart( : )  ) ;  end
  if  ~ isscalar (  mstop   )  ,  mstop  = max (  mstop ( : )  ) ;  end
  
  % If mstart signal and analysis epoch is missing from NSP recording , if
  % mstart signal missing and we're using header data to estimate peaks ,
  % or the mstop signal is missing
  if  ( isempty( mstart )  &&  ~ EPOCHF )  ||  ...
      ( isempty( mstart )  &&  par.useheader )  ||  ...
      ( isempty( mstop  )  && ~par.useheader )
    n2perr = EMISSING ;
    return
  end % check for missing data
  
  
  %%% Search for photodiode peaks %%%
  
  % Search for specified electrode
  i = [ nsx.ElectrodesInfo.ElectrodeID ]  ==  par.eid ;
  
  % Not found
  if  all (  ~ i  )
    
    error (  'MAK:maknsp2ptb:eid'  ,  [ 'maknsp2ptb: cannot find ' , ...
      'electrode id %d' ]  ,  par.eid  )
    
  end % no electrode
  
  % Convert photodiode measurements to double floating point numbers
  ph = double (  nsx.Data( i , : )  ) ;
  
  % Maxima
  [ p.max , t.max ] = findpeaks (  + ph ,  nsx.MetaTags.SamplingFreq ,  ...
    'MinPeakHeight' ,  0 ,  'MinPeakProminence' ,  par.minprom ,  ...
    'MinPeakWidth' ,  par.minwid  ) ;
  
  % Minima
  [ p.min , t.min ] = findpeaks (  - ph ,  nsx.MetaTags.SamplingFreq ,  ...
    'MinPeakHeight' ,  0 ,  'MinPeakProminence' ,  par.minprom ,  ...
    'MinPeakWidth' ,  par.minwid  ) ;
  
  % Restore minima
  p.min = - p.min ;
  
  % findpeaks assumes that the first sample is at time zero but we will add
  % one sample's duration
  t.max = t.max  +  1 / nsx.MetaTags.SamplingFreq ;
  t.min = t.min  +  1 / nsx.MetaTags.SamplingFreq ;
  
  
  %%% Sample tail end off-trial peaks %%%
  
  % Compute classification limits of maxima and minima in photodiode peaks
  % following the end of the trial
  for  C = { { 'min' , 'tailburnmin' } , { 'max' , 'tailburnmax' } }
    
    % Map to generic names , field name and tail burn
    [ f , tb ] = C{ 1 }{ : } ;
    
    % Use data before the trial starts
    if  par.useheader
      
      % Locate peaks before the start of the trial
      i = find (  mstart  >  t.( f )  ) ;
      
      % Discard burn-in peaks
      if  ~ isempty (  i  )
        i( end - par.( tb ) + 1 : end ) = [] ;
      end
      
    % Using data following the end of the trial
    else
      
      % Locate peaks that follow the end of the trial
      i = find (  mstop  <  t.( f )  ) ;
      
      % Discard burn-in peaks
      if  ~ isempty (  i  )
        i( 1 : par.( tb ) ) = [] ;
      end
    
    end % heads or tails?
    
    % No peaks left , we can't estimate classification limits , return
    % error
    if  numel ( i )  <  3
      n2perr = ENOREGRESS ;
      return
    end
    
    % Calculate classification limits
    clim.( f ) = median (  p.( f )( i )  )  +  ...
      [ - par.stdevs , + par.stdevs ] * std (  p.( f )( i )  ) ;
    
  end % find class limits
  
  
  %%% Look for head end off-trial peaks %%%
  
  % Use classification boundaries required for detecting the synchronising
  % flip
  f = par.sync ;
  
  % Find flips prior to the mstop signal that lie within classification
  % boundaries. If using header then don't bother about mstop signal.
  i = clim.( f )( 1 ) <= p.( f )  &  p.( f ) <= clim.( f )( 2 ) ;
  
  if  ~ par.useheader  ,  i = t.( f ) < mstop  &  i ;  end
  
  % There must be a chain of peaks at the head of the recording that lie
  % within the classification boundary. Find the last of them.
  i = i( 1 )  *  find (  diff( i )  ,  1  ,  'first'  ) ;
  
  % There are no off-trial peaks at the start of the trial and the analysis
  % epoch is missing
  if  isempty( i )  ||  i == 0
    
    % However, if the analysis epoch is recorded then continue
    if  EPOCHF
    
      % Guarantee i is zero in this case
      i = 0 ;
      
    else
      
      % Quit with error
      n2perr = EMISSING ;
      return
      
    end

  end % no starting off-trial peaks
  
  
  %%% Now detect the synchronising flip %%%
  
  % Make sure that the sync parameter is valid
  if  ~ any ( strcmp(  par.sync  ,  { 'min' , 'max' }  ) )
    
    error (  'MAK:maknsp2ptb:par'  ,  ...
      'maknsp2ptb:par.sync has invalid string ''%s'''  ,  par.sync  )
    
  end % sync parameter check
  
  % Determine which operator (fh) and classification boundary (cb) is
  % required to find the sync flip
  switch  par.sync
    case  'min'  ,  fh = @lt ;  cb = clim.( f )( 1 ) ;
    case  'max'  ,  fh = @gt ;  cb = clim.( f )( 2 );
  end
  
  % The header off-trial peaks will come just before the sync flip. Hence,
  % start searching from the last header off-trial peak. Only do this if i
  % is non-zero.
  if  0  <  i
    
    for  i = i + 1 : numel( p.( f ) )

      % See if the flip exceeds the off-trial classification boundary in
      % the desired direction. Break loop when it is found.
      if  fh (  p.( f )( i )  ,  cb  )  ,  break  ,  end

    end % sync flip
  
  end % is i non-zero?
  
  % If not found then quit on error
  if  i  ==  numel( p.( f ) )
    n2perr = ENOREGRESS ;
    return
  end
  
  % Get synchronising flip time if i is non-zero
  if  0  <  i
    
    tsync = t.( f )( i ) ;
    
  % Otherwise, zero will be the start time. We only get here if data from
  % the start of the recording is missing while the analysis epoch was
  % recorded and skipped frames are allowed. Cheat by 
  else
    
    % Zero start time
    tsync = 0 ;
    
  end
  
  
  %%% Detect in-trial flips %%%
  
  % Use specified classification boundaries
  f = par.flips ;
  
  % Return classification limits for indentifying in-trial flips
  n2pclim = clim.( f ) ;
  
  % Find all peaks that are outside of classification boundaries and come
  % after the sync flip
  i = ( p.( f ) < clim.( f )( 1 )  |  p.( f ) > clim.( f )( 2 ) )  &  ...
    tsync  <  t.( f ) ;
  
  % Return photometer measurements
  n2pvolts = p.( f )( i ) ;
  
  % Time each flip from peak that immediately preceeds those identified
  if  par.leadpeak
    
    % Determine which data set to use
    switch  f
      case  'min'  ,  f = 'max' ;
      case  'max'  ,  f = 'min' ;
    end
    
    % Get times of preceeding peaks
    tframe = precpeaks (  t.( par.flips )( i ) ,  t.( f )  ) ;
    
  else
    
    % Get measured frame onset times
    tframe = t.( f )( i ) ;
    
  end % leadpeak
  
  % No frame times were found , quit with error
  if  isempty (  tframe  )
    n2perr = ENOREGRESS ;
    return
  end
  
  % There are skipped frames , or start-of-trial data is missing but
  % analysis epoch is not (this is true when tsync == 0)
  if  par.skips  &&  (  any( ptb.Missed )  ||  tsync == 0  )
    
    % Look for number consecutive non-skipped frames , from head ...
    nsh = find (  ptb.Missed  ,  1  ,  'first'  )  -  1 ;
    
    % ... and tail
    nst = numel (  ptb.Missed  )  -  ...
      find (  ptb.Missed  ,  1  ,  'last'   ) ;
    
    % Missing data but analysis epoch present , adjust skip-frame numbers
    if  tsync  ==  0
      
      % Discard all header frames
      nsh = 0 ;
      
      % Keep as many tail frames as possible
      if  isempty( nst )
        nst = min ( [  numel( tframe )  ,  numel( ptb.Missed )  ] ) ;
      else
        nst = min ( [  numel( tframe )  ,  nst  ] ) ;
      end
      
    end % missing data
    
    % Not enough non-skipped frames
    if  nsh  +  nst  <  par.skipn
      
      n2perr = ENOREGRESS ;
      return
    
    end % not enough non-skipped
    
    % Get non-skipping measured and reported frame onset times
    tframe(  nsh + 1 : end - nst  ) = [] ;
    ptbson = ptb.StimulusOnsetTime( [  1 : nsh ,  end - nst + 1 : end  ] );
    
  % No skipped frames , if the number of measured frame onsets is not equal
  % to those reported by PTB then there is missing data somewhere. Error
  % and quit. Not a problem if relaxed matching enabled.
  elseif  ~ par.relaxedmatch  &&  ...
      numel (  tframe  )  ~=  numel (  ptb.StimulusOnsetTime  )
    
    n2perr = ENOREGRESS ;
    return
    
  % No skipped frames and the number of measured onsets equals the number
  % of reported ones , point to PTB frame onsets
  else
    
    % Point to PTB record of stimulus onset times
    ptbson = ptb.StimulusOnsetTime ;
    
    % Look for the minimum number of frames
    mi = min( numel( tframe ) , numel( ptbson ) ) ;
    
    % Eliminate excess frames from one set or the other
    ptbson( mi + 1 : end ) = [ ] ;
    tframe( mi + 1 : end ) = [ ] ;
    
  end % skips
  
  % Convert PTB times to seconds in double floating point , from number of
  % microseconds in integers
  ptbson = double (  ptbson  )  /  1e6 ;
  
  
  %%% Compute regression %%%
  
  n2pcoef = robustfit (  tframe( : )  ,  ptbson( : )  )' ;
  
  
end % maknsp2ptb


%%% Subroutines %%%

% Find times from set t2 that preceeds each value in set t1
function  t = precpeaks (  t1  ,  t2  )
  
  % Build one grand vector of data , with accompanying identity vector
   t = [  t1  ,  t2  ] ;
  id = [  ones( size( t1 ) )  ,  2 * ones( size( t2 ) )  ] ;
  
  % Sort time stamps , ascending
  [ t , j ] = sort (  t  ) ;
  
  % Apply sorting to id vector
  id = id( j ) ;
  
  % Return time stamps that belong to set 2 and come just before each time
  % from set 1
  t = t (  id( 1 : end - 1 ) == 2  &  id( 2 : end ) == 1  ) ;
  
end % precpeaks

