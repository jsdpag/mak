
function  [ outcome , msig , minlatflg ] = ...
   makreclasstrial (  p ,  MCC ,  MSID ,  sd ,  td ,  msig ,  eye ,  hit  )
% 
% outcome = makreclasstrial( p , MCC , MSID , sd , td , msig , eye , hit )
% 
% MET Analysis Kit, pre-processing. A trial's outcome may have been mis-
% classified as broken in cases where the binocular eye signal failed to
% enter a hit region, although the subject was looking at the associated
% stimulus. These trials will be re-classified by allowing either monocular
% eye signal to indicate the target. The task logic will be used to
% determine the actual outcome. A trial outcome character is returned in
% outcome: 'b' for broken, 'c' for correct, and 'f' for failed. If the
% trial is reclassified then an updated version of msig is returned in
% which the mtarget signal is assigned a new cargo value corresponding to
% the monocularly targeted stimulus.
% 
% Broken trials are candidates for reclassification if the mtarget signal
% leading to the break occurred at least some minimum time after the onset
% of a reference state in the task logic; if this occurred then minlatflg
% returns true, otherwise it returns false to say that the minimum latency
% was not reached. A trial is then reclassified if only one of its
% monocular eye signals fell within a hit region at the time of the mtarget
% signal. Otherwise, the trial remains broken. If the reference state is
% not present in the task logic used by the trial then 'b' is returned.
% 
% 
% Input arguments
% 
% p - Is a struct with at least the fields .state and .minlatency. The
%   first names a task logic state. The second gives a latency relative to
%   the onset of that state, in seconds.
% 
% MCC - MET controller constants, as returned by metctrlconst
% 
% MSID - Struct mapping MET signal names to signal identifiers.
% 
% sd - Session descriptor loaded from the session's sessdesc.mat file.
% 
% td - Trial descriptor loaded from the trial's param_*.mat file.
% 
% msig - Struct containing MET signal variables loaded from the trial's
%   metsigs_*.mat file. 
% 
% eye - Struct containing eye positions. The 'eye' variable loaded from
%   the trial's eyepos_*.mat file.
% 
% hit - Struct containing hit region variables loaded from the trial's
%   hitregion_*.mat file.
% 
% 
% NOTE: Tests whether the reference state has timed out by comparing the
%   break time against onset of the state. This is wrong. The correct way
%   is to examine the record of PTB frame times to find the expected onset
%   of the next frame at the time that the mtarget signal was received.
%   However, the simple method should work most of the time, while good
%   task logic design will make the problem irrelevant i.e. the same end
%   state should be reached from the reference state before and after
%   timeout for a given targeted task stimulus.
%
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Constants %%%
  
  % Column indices for eye position coordinates
  C = struct (  'LEFT_XY' ,  [ 1 , 2 ]  ,  'RIGHT_XY' ,  [ 3 , 4 ]  ) ;
  
  
  %%% Default return value %%%
  
  % In the event that the trial really is broken, return 'b' by default
  outcome = 'b' ;
  
  % Minimum latency not reached
  minlatflg = false ;
  
  
  %%% Find time of minimum latency %%%
  
  % Get the task logic used by this trial
  l = sd.logic.( td.logic ) ;
  
  % Reference state is not in the task logic. Quit here.
  if  all ( ~ strcmp(  p.state  ,  l.nstate  ) )  ,  return  ,  end
  
	% Get reference state's identifier
  istate = l.istate.( p.state ) ;
  
  % Return the onset time of the final occurrence of that state
  reftim = sigtim (  msig  ,  MSID.mstate  ,  istate  ) ;
  
  % Add the relative latency to get the absolute time of earliest allowable
  % break
  minlat = reftim  +  p.minlatency ;
  
  % State did not occur , quit
  if  isempty (  minlat  )  ,  return  ,  end
  
  
  %%% Find time of break %%%
  
  % Get time of the final mtarget signal reporting that the subject
  % targeted 'none' i.e. targeted no stimulus on screen. And get its 
  [ tnone , imtsig ] = sigtim (  msig  ,  MSID.mtarget  ,  l.istim.none  );
  
  % No mtarget/none signal was found , this happens if the trial times out.
  % Or the break occurred prior to the minimum latency. Quit either way.
  if  isempty( tnone )  ||  tnone < minlat  ,  return  ,  end
  
  % Minimum latency has been reached
  minlatflg = true ;
  
  
  %%% Get eye position at time of break %%%
  
  % Find the eye sample that occurred as close as possible to the breaking
  % mtarget signal
  [ ~ , i ] = min ( abs(  eye.time  -  tnone  ) ) ;
  
  % Convert eye positions to single floating point values in units of
  % degrees of visual field
  pos = double (  eye.position( i , : )'  )  /  100 ;
  
  
  %%% Hit regions at time of break %%%
  
  % Find hit region times that are up to the break
  i = hit.time  <=  tnone ;
  
  % Get hit region reports up to that time
  H = hit.hitregion( i , : ) ;
  
  % Return the final set of hit regions for each stimulus link
  H = cellfun (  @getlasthit  ,  num2cell( H , 1 )  ,  ...
    'UniformOutput'  ,  false  ) ;
  
  
  %%% Compare monocular positions to hit regions %%%
  
  % Left eye
  lh = cellfun(  @( h ) makgethits( MCC , h , pos( C.LEFT_XY  ) )  ,  H  );
  
  % Right eye
  rh = cellfun(  @( h ) makgethits( MCC , h , pos( C.RIGHT_XY ) )  ,  H  );
  
  % Exclusive OR to find hit regions targeted by only one eye
  h = xor (  lh  ,  rh  )  ;
  
  % Quit now when both monocular eye positions are outside of all hit
  % regions , or when they are in separate hit regions
  if  sum (  h  )  ~=  1  ,  return  ,  end
  
  
  %%% Reclassify the trial %%%
  
  % Name of task used by this trial
  ntask = td.task ;
  
  % Name of task's stimulus link that was targeted monocularly
  nlink = hit.stimlink{ h } ;
  
  % Get name of task-logic stimulus that was targeted monocularly.
  % Remember, the task-logic stimulus name is linked via a task's stimulus
  % link to a MET ptb stimulus definition. So we had to work backwards from
  % definition to task-logic name.
  nstim = sd.task.( ntask ).link.( nlink ).stim ;
  
  % The task logic identifier of the named task-logic stimulus
  istim = l.istim.( nstim ) ;
  
  % Determine whether the state had timed out (2) or not (1). This is the
  % simple way, but prone to error. To do this properly requires looking at
  % the PTB frame timestamp record. For now, we will assume that this is
  % correct. Good task logic design plus the infrequency of this being
  % wrong should guard against incorrect reclassification.
  tout = ( l.T( istate )  <=  tnone - reftim )  +  1 ;
  
  % Identifier of the end state
  i = l.E(  istate  ,  istim  ,  tout  ) ;
  
  % Not an end state? This can happen for at least two reasons. First, we
  % finally mis-judged if the state had timed out and used the wrong lookup
  % table. Second, the subject's eyes wandered off of a legitimate target
  % and into the void. In either case we can recover the pre-processing
  % session by simply rejecting this trial.
  if  i  <=  l.N.state - 4  ,  return  ,  end
  
  % Reclassify the outcome of the trial based on the name of the end state
  switch  l.nstate{ i }
    case  'correct'  ,  outcome = 'c' ;
    case  'failed'   ,  outcome = 'f' ;
  end % reclass
  
  % Update cargo value of mtarget signal
  msig.crg( imtsig ) = istim ;
  
  
end % makreclasstrial


%%% Sub-routines %%%

% Returns the time of the final MET signal to have signal identifier sid
% and cargo value crg
function  [ t , i ] = sigtim (  msig  ,  sid  ,  crg  )
  
  % Look for all occurrences of this state
  i = find (  msig.sig == sid  &  msig.crg == crg  ) ;
  
  % Get maximum onset time and signal index
  [ t , j ] = max ( msig.tim(  i  ) ) ;
  
  % Return signal index
  i = i (  j  ) ;
  
end % sigtim


% Returns the last non-empty element of H
function  h = getlasthit ( H )
  
  % Find the last non-empty element
  i = find (  ~ cellfun( @isempty , H )  ,  1  ,  'last'  );
  
  % Return it
  h = H{ i } ;
  
end % getlasthit

