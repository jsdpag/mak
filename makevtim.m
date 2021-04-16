
function  t = makevtim (  d  ,  typ  ,  val  )
% 
% t = makevtim (  d  ,  typ  ,  val  )
% 
% 
% MET Analysis Kit. Returns vector t of event times relative to the start
% of the analysis epoch for events of type typ with value val for all
% trials in makprep data struct d. The last instance of an event in each
% trial is used.
% 
% typ is a single character code saying what type of event to search for:
% 
%   s - Change of task logic state
%   t - New target selected by subject
%   r - Reward given
%   y - Reward type
% 
% val is a valid value for the event data, dependent on the type of event:
%   
%   typ:
%   
%   's' - val must be a string naming the task logic state
%   't' - val must be a string naming the task stimulus
%   'r' - val must be a single numeric value giving the duration in
%     milliseconds of the reward
%   'y' - val must be a non-zero natural number identifying the reward type
%     e.g. 1, 2, 3, etc.
%   
%   See the task logic for state and stimulus names. d.sd( 1 ).logic
%   contains a description of all task logics used in the experiment.
%   d.sd( 1 ).( lname ).nstate and d.sd( 1 ).( lname ).nstim give the state
%   and stimulus names for logic lname, respectively.
% 
% val can be empty i.e. [ ] so that the final event of specified type is
% used, regardless of the event data.
% 
% t is a double row vector with d.numtrials elements. If no event can be
% found for trial i then t( i ) is NaN.
% 
% 
% Written by Jackson Smith - January 2019 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % Fields required by d
  DFIELDS = {  'sd'  ,  'tasknames'  ,  'logicnames'  ,  'task_ind'  ,  ...
    'logic_ind'  ,  'numtrials'  ,  'event'  } ;
  
  % Fields required by d.event
  EFIELDS = {  'time'  ,  'type'  ,  'data'  } ;
  
  % Valid characters for typ
  TYPCHAR = {  's'  ,  't'  ,  'r'  ,  'y'  } ;
  
  
  %%% Input check %%%
  
  % Number of inputs
   narginchk (  3  ,  3  )
  nargoutchk (  0  ,  1  )
  
  % Is d a single element struct?
  if  ~ isstruct (  d  )  ||  ~ isscalar (  d  )
    
    error (  'MAK:makevtim:dstruct'  ,  ...
      'makevtim: d must be a single element struct'  )
  
  % Does d have the required fields?
  elseif  ~ all ( isfield(  d  ,  DFIELDS  ) )
    
    error (  'MAK:makevtim:dfields'  ,  ...
      'makevtim: d must have fields: %s'  ,  strjoin( DFIELDS , ' , ' )  )
    
  % Check that event is a single element struct
  elseif  ~ isstruct (  d.event  )  ||  ~ isscalar (  d.event  )
    
    error (  'MAK:makevtim:estruct'  ,  ...
      'makevtim: d.event must be a single element struct'  )
    
  % Does d.event have required fields?
  elseif  ~ all ( isfield(  d.event  ,  EFIELDS  ) )
      
    error (  'MAK:makevtim:efields'  ,  ...
      'makevtim: d.event must have fields: %s'  ,  ...
        strjoin( EFIELDS , ' , ' )  )
      
	% typ must be one of these characters
  elseif  ~ any ( strcmp(  typ  ,  TYPCHAR  ) )
    
    error (  'MAK:makevtim:typ'  ,  ...
      'makevtim: typ must be one of chars: %s'  ,  ...
        strjoin( TYPCHAR , ' , ' )  )
    
  end % input check
  
  % Note - at this point we have not checked the number of trials, or that
  % vectors of fields in d or d.event have a matching number of elements.
  % However, if d has passed the tests presented then we will assume it was
  % produced by makprep and possibly loaded by makload.
  
  % Don't check val if it is empty
  if  ~ isempty (  val  )
    
    % Checking val depends on the value of typ , unless val is empty
    switch  typ

      % State or stimulus name
      case  { 's' , 't' }

        % Must be a string
        if  ~ ischar (  val  )  ||  ~ iscellstr (  { val }  )

          error (  'MAK:makevtim:val'  ,  ...
            'makevtim: for typ ''s'' or ''t'' , val must be a string'  )

        end % check string

      % Reward duration
      case  'r'

        % May be empty
        if  isempty (  val  )

          % Does nothing other than to skip the next test

        % Scalar non-negative integer value
        elseif  ~ isscalar (  val  )  ||  ~ isnumeric (  val  )  ||  ...
            val  <  0  ||  mod (  val  ,  1  )

          error (  'MAK:makevtim:val'  ,  ...
            'makevtim: for typ ''r'' , val must be a non-negative integer'  )

        end % reward duration

      % Reward type
      case  'y'

        % Scalar non-zero natural number
        if  ~ isscalar (  val  )  ||  ~ isnumeric (  val  )  ||  ...
            val  <=  0  ||  mod (  val  ,  1  )

          error (  'MAK:makevtim:val'  ,  [ 'makevtim: for typ ''y'' , ' ,...
            'val must be a non-zero natural number' ]  )

        end % reward duration

      % Unrecognised type
      otherwise

        error (  'MAK:makevtim:typ_unrecognised'  ,  ...
        'makevtim: typ unrecognised in val input check'  )

    end % check val
  
  end % val not empty
  
  
  %%% Event times %%%
  
  % Flag saying whether we need to check event data , no need if val is [ ]
  flgval = ~ isempty (  val  ) ;
  
  % Task logic field used to get either the state or stimulus identifier.
  % Only needed if we're checking event data value.
  if  flgval
    
    switch  typ
      case  's'  ,  sidfld = 'istate' ;
      case  't'  ,  sidfld = 'istim'  ;
      otherwise  ,  sidfld = '' ;
    end % state/stim identifier field
  
  end % chk ev dat
  
  % Flag saying that we need to get the state or stimulus logic identifier.
  % This is automatically false if we are ignoring the event data value.
  flgsid = flgval  &&  ~ isempty (  sidfld  ) ;
  
  % Allocate output and initialise to NaN , the value returned if an event
  % is not found
  t = nan (  1  ,  d.numtrials  ) ;
  
  % Trials
  for  i = 1 : d.numtrials
    
    % Get the state or stimulus identifier as the search value
    if  flgsid
      
      % Task logic name
      ln = d.logicnames{  d.logic_ind( i )  } ;
      
      % Check that named state or stimulus exists in task logic
      if  ~ isfield (  d.sd( 1 ).logic.( ln ).( sidfld )  ,  val  )
        
        % Doesn't exist so we won't find it , t( i ) is NaN so skip to next
        % trial
        continue
        
      end % check field exists
      
      % The identifier
      v = d.sd( 1 ).logic.( ln ).( sidfld ).( val ) ;
      
    % Use numeric value given in val if we're checking event values
    elseif  flgval
      
      v = val ;
      
    end % search value
    
    % Event times, types, and data
    etim = d.event.time{ i } ;
    etyp = d.event.type{ i } ;
    
    % Find event of the specified type
    e = etyp  ==  typ ;
    
    % Check event data
    if  flgval
      
      % Event data
      edat = d.event.data{ i } ;
      
      % Find events with specified data value
      e = e  &  edat'  ==  v ;
      
    end % chk edat
    
    % Find last instance of matching event
    j = find (  e  ,  1  ,  'last'  ) ;
    
    % None found , t( i ) is NaN so skip to next trial
    if  isempty (  j  )  ,  continue  ,  end
    
    % Get event time
    t( i ) = etim( j ) ;
    
  end % trials
  
  
end % makevtim

