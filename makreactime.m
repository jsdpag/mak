
function  ...
[ rt , start , duration , angle , amplitude , x_start , y_start ] = ...
makreactime(  MCC ,  MSID ,  BC ,  par ,  sd ,  td ,  msig ,  eye ,  hit  )
% 
% [ rt , start , duration , angle , amplitude , x_start , y_start ] = 
%        makreactime ( MCC , MSID , BC , par , sd , td , msig , eye , hit )
% 
% MET Analysis Kit, pre-processing. Computes trial's reaction time.
%
% 
% Input arguments
% 
% BC - Struct with fields .b and .a containing the coefficients for a
%   lowpass butterworth filter.
% 
% MSID - Struct mapping MET signal names to signal identifiers.
% 
% par - Struct with reaction time parameters. Must have at least fields:
% 
%   .velocity - Velocity threshold in degrees per second
%   .accel - Acceleration threshold in degrees per second-squared
%   .coef - Multiplicative coefficient applied to thresholds
%   .state - Name of task-logic reference state, if empty or not found then
%      mstart time is used
%   .offsets - Two element vector giving start and end time of a window
%     relative to the reference state, for sampling velocity and
%     acceleration.
%   .lognvel - Area under log-normal curve fit to velocities to get
%     threshold.
%   .lognacc - Same but for acceleration.
%   .stim - Name of task-logic reference stimulus , ignored if empty
%   .verifythr - If non-zero then this is taken to be the number of seconds
%     prior to the final mtarget event in which to compare the velocity and
%     acceleration values against the measured thresholds.
%   .halfmax - Logical flag. If non-zero then .verifythr defines an
%     analysis window prior to the final mtarget event. This window is
%     searched with the .logn* parameters to define significantly high
%     velocity and acceleration samples. The peak of these significant
%     samples is found. Half the peak values are then used as the
%     thresholds. What if the .logn*/.verifythr thresholds fail to identify
%     significant values? Then search the .verifythr window for the maximum
%     value and return half of that.
% 
% sd , td - Session and trial descriptors
% 
% msig - Struct with MET signal variables loaded from trial's metsigs_*.mat
%   file 
% 
% eye - Eye position samples and times from 'eye' variable loaded from
%   trial's eyepos_*.mat file 
% 
% hit - Struct with hit position variables loaded from trial's
%   hitregion_*.mat file
% 
% 
% Output arguments , all are scalar , default zero returned if unable to
% calculate value
% 
% rt - Reaction time relative to reference state , in seconds
% 
% start - Time from start of the trial when saccade began , in seconds
% 
% duration - Duration of saccade , in seconds
% 
% angle - Polar coordinate angle of saccade , in degrees
% 
% amplitude - Polar coordinate radius (saccad amplitude) , in degrees of
%   visual field
% 
% x_start - Eye position x-coordinate at the start of the saccade
% 
% y_start - Eye position y-coordinate at the end of the saccade
% 
% 
% Reaction time is defined as the time that the last saccade onto a
% stimulus began, relative to the onset of a reference state. This is
% computed by first averaging together binocular eye positions, then
% filtering them with a low-pass butterworth filter. Velocity and
% acceleration are computed from the averaged, filtered eye positions.
% The velocity and acceleration thresholds are used to determine when the
% eyes were fixating. Fixations may occur sacmin seconds before the final
% mtarget signal, at latest. Another constraint is that the eyes must be
% looking at a specified reference stimulus at the time of pre-saccadic
% fixation, such as a fixation target. Lastly, the start of the saccade
% must follow at least fixdur seconds of fixating eye positions that are
% below both thresholds. If the start of a saccade cannot be found, for
% instance if there is too much noise in the eye signal, then the time
% of the final eye position within the hit region of the reference
% stimulus is used ; the final mtarget time is used if no reference
% stimulus is provided or if the averaged-filtered eye position never left
% the reference stimulus.
% 
% If a reference state is provided, then an analysis window relative to
% this is used to sample eye velocity and acceleration. Log-normal curves
% are fit to each distribution and used to estimate the velocity and
% acceleration thresholds. Areas under the curve are given in lognvel and
% lognacc so that the inverse of the fitted curves provides appropriate
% thresholds. When no reference state is provided then the manual values
% are used.
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Constants %%%
  
  % Left and right eye position column indices for x-axis
  X = [ 1 , 3 ] ;
  
  % Left and right eye position column indices for y-axis
  Y = [ 2 , 4 ] ;
  
  % Eye sampling rate
  EYESHZ = MCC.SHM.EYE.SHZ ;
  
  % Set thresholds
  SVTH = par.coef  *  par.velocity ;
  SATH = par.coef  *  par.accel    ;
  
  % Use these thresholds by default
  VTH = SVTH ;
  ATH = SATH ;
  
  % Minimum number of fixation samples required
  FIXDUR = ceil (  par.fixdur  *  EYESHZ  ) ;
  
  
  %%% Default return values %%%
  
  start = 0 ;
  duration = 0 ;
  angle = 0 ;
  amplitude = 0 ;
  x_start = 0 ;
  y_start = 0 ;
  
  
  %%% Trial event times %%%
  
  % Locate mstart time
  mstart = msig.tim(  msig.sig  ==  MSID.mstart  ) ;
  
  % No mstart time , raise error
  if  isempty ( mstart )

    error (  'MAK:makreactime:nomstart'  ,  [ 'makreactime: ' , ...
      'no mstart signal for %s, exp %d, sess %d, trial %d' ]  ,  ...
      sd.subject_id  ,  sd.experiment_id  ,  sd.session_id  ,  ...
      td.trial_id  )

  end % no mstart
  
  % Last mtarget signal
  i = find (  msig.sig  ==  MSID.mtarget  ,  1  ,  'last'  ) ;
  mtarg = msig.tim( i ) ;
  
  % No mtarget found
  if  isempty ( mtarg )
    
    % Try mstop time
    i = find (  msig.sig  ==  MSID.mstop  ,  1  ,  'last'  ) ;
    mtarg = msig.tim( i ) ;
    
    % No mstop time , raise error
    if  isempty ( mtarg )
      
      error (  'MAK:makreactime:nomstart'  ,  [ 'makreactime: ' , ...
        'no mstop signal for %s, exp %d, sess %d, trial %d' ]  ,  ...
        sd.subject_id  ,  sd.experiment_id  ,  sd.session_id  ,  ...
        td.trial_id  )
      
    end % no mstart
    
  end % mstop
  
  % Reference state time initialised empty
  reftime = [] ;
  
  % Look for reference state if it is in trial's task logic
  if  ~ isempty (  par.state  )  &&  ...
      any ( strcmp(  par.state  ,  sd.logic.( td.logic ).nstate  ) )
    
    % Reference state identifer
    i = sd.logic.( td.logic ).istate.( par.state ) ;
    
    % Search for last onset of reference state
    i = find (  msig.sig == MSID.mstate  &  msig.crg == i  ,  ...
      1  ,  'last'  ) ;
    
    % Get state onset time if we found it
    if  ~ isempty ( i )  ,  reftime = msig.tim( i ) ;  end
    
  end % ref state
  
  % Reference state time is still empty
  if  isempty (  reftime  )
    
    % Empty the reference state parameter to signal that manual thresholds
    % should be used. Copy-on-write will protect the original value in the
    % calling function.
    par.state = '' ;
    
    % Use mstart time instead
    reftime = mstart ;
    
  end % use mstart
  
  
  %%% Reference fixation stimulus %%%
  
  % Reference stim hit regions , initialise empty
  hrefstim = [] ;
  
  % Reference stimulus provided and it is in task logic
  if  ~ isempty (  par.stim  )  &&  ...
      any ( strcmp(  par.stim  ,  sd.logic.( td.logic ).nstim  ) )
    
    % Look for stimulus links onto this task stimulus
    i = strcmp (  par.stim  ,  { td.stimlink.nstim }  ) ;
    
    % Get corresponding hit regions
    hrefstim = hit.hitregion( : , i ) ;
    
    % Find hit regions up to the last mtarget time
    i = hit.time  >  mtarg ;
    hrefstim( i ) = [] ;
    
    % Take the last hit region provided for each link. Here we split each
    % column into its own cell array. The index of the final non-empty
    % element is found and used to retrieve the same element. All retrieved
    % elements are stored in a new cell array on return.
    hrefstim = cellfun (  ...
      @( h ) h{ find( ~cellfun( @isempty , h ) , 1 , 'last' ) }  ,  ...
        num2cell( hrefstim , 1 )  ,  'UniformOutput'  ,  false  ) ;
    
  end % reference stim
  
  
  %%% Eye position signal processing %%%
  
  % Turn eye positions into double precision floating point numbers , then
  % convert from centi-degrees to degrees
  p = double (  eye.position  )  /  100 ;
  
  % Average binocular positions , column order is [ x-axis , y-axis ]
  % coordinates
  p = [  mean(  p( : , X )  ,  2  )  ,  mean(  p( : , Y )  ,  2  )  ] ;
  
  % Apply low-pass filter , if requested
  if  par.lowpass  ,  p = filtfilt ( BC.b , BC.a , p ) ;  end
  
  % Compute velocity marginals ...
  vx = diff( p( : , 1 ) )  *  EYESHZ ;
  vy = diff( p( : , 2 ) )  *  EYESHZ ;
  
  % ... and vector , zero padded at head
  v = [  0  ;  sqrt(  vx .^ 2  +  vy .^ 2  )  ] ;
  
  % Compute acceleration marginals ...
  ax = diff ( vx )  *  EYESHZ ;
  ay = diff ( vy )  *  EYESHZ ;
  
  % ... and vector , zero padded at head. Zero padding keeps data in
  % register with eye positions and sample times.
  a = [  0  ;  0  ;  sqrt(  ax .^ 2  +  ay .^ 2  )  ] ;
  
  % Find samples within reference stimulus hit region
  if  ~ isempty (  hrefstim  )
    
    % Returns cell vector of logical row vectors
    irefstim = cellfun (  @( H ) makgethits( MCC , H , p' )  ,  ...
      hrefstim  ,  'UniformOutput'  ,  false  ) ;
    
    % Element-wise OR into a single logical row vector and transpose
    irefstim = any (  cell2mat( irefstim )  ,  1  )' ;
    
    % Look for eye positions from the reference time to the start of the
    % shortest allowable saccade ...
    i = reftime <= eye.time  &  eye.time < mtarg ;
    
    % ... and examine the case when the average, filtered eye position
    % never leaves the reference stimulus hit region. This could happen if
    % the gaze location of only one eye falls off the edge. When this
    % happens, return the mtarget time minus the reference time as the
    % reaction time ... the best we can do in this case.
    if  all (  irefstim(  i  )  )
      rt = mtarg  -  reftime ;
      return
    end
    
  % No reference stimulus , use all eye positions
  else
    
    irefstim = true (  size(  p  ,  1  )  ,  1  ) ;
    
  end % ref stim vs eyes
  
  
  %%% Measure thresholds %%%
  
  % If reference state provided and we're allowed to measure a threshold
  if  ~ isempty (  par.state  )  &&  ...
      (  0  <=  par.lognvel  ||  0  <=  par.lognacc  )
    
    % Find eye samples within analysis time window relative to a reference
    % state
    i = reftime + par.offsets( 1 ) <= eye.time  &  ...
        reftime + par.offsets( 2 ) >= eye.time ;

    % Refine this. Drop three samples prior to the eyes leaving the
    % fixation target and after returning. This is to discard some
    % samples when the eyes were still accelerating. Altogether , this
    % should help to ignore blink artefacts ... but see verification step
    % below.
    i = chopaccel (  i  ,  irefstim  ) ;
    
    % Measure velocity threshold
    if  0  <=  par.lognvel
      
      VTH = measureth (  par.lognvel  ,  v( i )  ) ;
    
    end % measure vel thresh
    
    % Same again for acceleration
    if  0  <=  par.lognacc
      
      ATH = measureth (  par.lognacc  ,  a( i )  ) ;
      
    end % measure accel thresh
    
    % Sanity checks. It is possible that blinks caused big velocity and
    % acceleration artefacts during the fixation period that were ignored
    % by MET because of the blink filter. But this might cause our
    % estimated velocity and acceleration thresholds to become very large.
    % If threshold verification is enabled then check current measured
    % thresholds before estimating new ones.
    if  par.verifythr
      
      % Find eye samples within a certain time of the mtarg event
      j = mtarg - par.verifythr  <=  eye.time  &  mtarg  >=  eye.time ;
      
      % Check velocity threshold
      if  0  <=  par.lognvel
        
        VTH = verifythr (  VTH  ,  v  ,  i  ,  j  ,  SVTH  ) ;
      
      end % vel thr
      
      % Check acceleration threshold
      if  0  <=  par.lognacc
        
        ATH = verifythr (  ATH  ,  a  ,  i  ,  j  ,  SATH  ) ;
      
      end % vel thr
      
      % Half-maximum flag is up. Use the thresholds we've got to identify
      % significantly high velocity and acceleration samples. Locate the
      % maximum of each, and return the half-maximum as the threshold we
      % want.
      if  par.halfmax
        
        % Half-max velocity threshold
        VTH = halfmaxthr( VTH , v( j ) ) ;
        
        % Half-max acceleration threshold
        ATH = halfmaxthr( ATH , a( j ) ) ;
        
      end % half-max thresholds
      
    end % verify thresholds
    
    
  end % measure thresholds
  
  
  %%% Search for end of saccade %%%
  
  % Look backwards from the mtarget time for the end of the saccade. This
  % must be when either velocity or acceleration is over threshold but
  % decreasing.
  i = eye.time < mtarg  &  ( ...
      (  v > VTH  &  diff( [ v ; 0 ] ) < 0  )  |  ...
      (  a > ATH  &  diff( [ a ; 0 ] ) < 0  )  ) ;
  
  % Look for the last sample when this is true
  j = find (  i  ,  1  ,  'last'  ) ;
  
  % End of saccade found
  if  ~ isempty (  j  )
    
    % Re-adjust par.sacmin to start saccade search from end of saccade
    par.sacmin = mtarg  -  eye.time( j ) ;
    
  end % end of saccade
  
  % If end of saccade not detected then default sacmin value is used
  
  
  %%% Search for start of saccade %%%
  
  % Find eye samples that preceed mtarget, are below thresholds, and on the
  % reference stimulus. In other words, these are candidate fixation
  % samples.
  i = eye.time < mtarg - par.sacmin  &  v < VTH  &  a < ATH  &  irefstim  ;
  
  % Initialise index of last fixation before saccade
  j = numel (  i  ) ;
  
  % Search for the last stretch of fixation prior to saccade
  while  j
    
    % Get the last candidate fixation
    j = find (  i( 1 : j )  ,  1  ,  'last'  ) ;
    
    % Nothing found , break loop
    if  isempty (  j  )  ,  j = 0 ;  break  ,  end
    
    % Grab a snippet of the eye classifications ending at the last
    % candidate fixation
    k = i (  max( [  1  ,  j - FIXDUR + 1  ] )  :  j  ) ;
    
    % All samples are fixations , we found the start of the saccade
    if  all (  k  )  ,  break  ,  end
    
    % There are possible saccadic samples in the snippet. Find the earliest
    % one and jump back to the sample just before it.
    j = j  -  FIXDUR  +  find (  ~ k  ,  1  )  -  1 ;
    
    % Make sure we don't have any negative index
    j = max ( [  0  ,  j  ] ) ;
    
  end % fix search
  
  % No saccade was found , return reaction time based on mtarget and return
  if  ~ j
    rt = mtarg  -  reftime ;
    return
  end
  
  
  %%% Measure saccade %%%
  
  % Store index of the first eye sample of saccade
  i = j  +  1 ;
  
  % First , get the onset time
  start = eye.time( i ) ;
  
  % Now compute reaction time
  rt = start  -  reftime ;
  
  % Get the starting point of the saccade
  x_start = p( j , 1 ) ;
  y_start = p( j , 2 ) ;
  
  % Find eye samples from the start of the saccade to the mtarget signal
  % that exceeds either threshold
  i = start <= eye.time  &  eye.time <= mtarg  &  ...
    ( VTH <= v  |  ATH <= a ) ;
  
  % And take the last one , plus one to find the first fixation sample
  % following saccade
  i = find (  i  ,  1  ,  'last'  )  +  1 ;
  
    % Don't allow i to exceed the last element of eye.time
    i = min (  i  ,  numel( eye.time )  ) ;
  
  % Duration of the saccade
  duration = eye.time( i )  -  start ;
  
  % Start of saccade from start of the trial
  start = start  -  mstart ;
  
  % Saccade end point
  x_end = p( i , 1 ) ;
  y_end = p( i , 2 ) ;
  
  % Compute angle
  opposite = y_end  -  y_start ;
  adjacent = x_end  -  x_start ;
  angle = atand (  opposite  /  adjacent  ) ;
  
    % If the saccade went left-wards then we must adjust the angle by
    % adding 180 degrees
    if  adjacent  <  0  ,  angle = angle  +  180 ;  end
    
    % Normalise angle between 0 and 360
    angle = mod (  angle  ,  360  ) ;
  
  % Amplitude is the Euclidean distance
  amplitude = sqrt( sum( [  adjacent  ,  opposite  ] .^ 2 ) ) ;
  
  
end % makreactime


%%% Sub-routines %%%

% Find eye samples within a time window and within range of a fixation
% stimulus. Drop 3 samples before the eyes leave the stimulus and 3 samples
% after they return. Why 3? Because it takes 3 samples to compute
% acceleration.
function  i = chopaccel (  t  ,  irefstim  )
  
  % Number of samples
  n = numel (  t  ) ;
  
  % Initialise output
  i = t  &  irefstim ;
  
  % Get indices of first and last eye sample within the time window
  f = find (  t  ,  1  ,  'first'  ) ;
  l = find (  t  ,  1  ,   'last'  ) ;
  
  % Find samples when eyes leave or return to fixation stimulus
  d = diff (  irefstim( f : l )  ) ;
  
  % Locate all samples when the eyes first leave the fixation stimulus
  j = f  +  find (  [ 0 ; d ]  ==  -1  )  -  1 ;
  
  % Remove 3 samples before each one
  for  j = j'  ,  i( max( [ 1 , j - 3 ] ) : j - 1 ) = 0 ;  end
  
  % Now locate all samples just before the eyes return to the fixation
  % stimulus
  j = f  +  find (  [ d ; 0 ]  ==  +1  )  -  1 ;
  
  % Remove 3 samples after each one
  for  j = j'  ,  i(  j + 1 : min( [ j + 3 , n ] )  ) = 0 ;  end
  
end % chopaccel


% Measure threshold from data x by estimating log-normal parameters , then
% taking the inverse value of the CDF at probability p
function  TH = measureth (  p  ,  x  )
  
  % Estimate mu and sigma parameters for log-normal distribution
     mu = mean ( log(  x  ) ) ;
  sigma =  std ( log(  x  ) ) ;

  % Estimate threshold
  TH = logninv (  p  ,  mu  ,  sigma  ) ;
  
end % measureth


% Verify measured input threshold TH against eye samples in x. If input TH
% exceeds the peak in x( j ) then estimate a new threshold from x( i )
% as the median plus the standard deviation. If that still exceeds the peak
% then return the set threshold STH.
function  TH = verifythr (  TH  ,  x  ,  i  ,  j  ,  STH  )
  
  % Find maximum value in x
  m = max (  x( j )  ) ;
  
  % Input threshold does not exceed the peak in x , so return input thresh
  if  TH  <=  m  ,  return  ,  end
  
  % Estimate a new threshold from the data
  TH = median (  x( i )  )  +  std (  x( i )  ) ;
  
  % New threshold does not exceed the peak in x , so return that
  if  TH  <=  m  ,  return  ,  end
  
  % Checks fail , return set threshold
  TH = STH ;
  
end % verifythr


% Use input TH threshold to define significant values of x. Max x is found
% and half the maximum is returned. If input TH fails to find significant
% values, then we skip that step and simply return half the max of x.
function  TH = halfmaxthr( TH , x )
  
  % Find significant data points
  i = TH  <=  x ;
  
  % We found some! Keep only significant data points.
  if  any( i )  ,  x = x( i ) ;  end
  
  % Return half the maximum
  TH = max( x )  ./  2 ;
  
end % halfmaxthr

