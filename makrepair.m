
function  makrepair( funstr , varargin )
% 
% makrepair( <function string> , <func arg 1> , <func arg 2> , ... )
% 
% MET Analysis Kit. Although one strives to obtain clean and well ordered
% data, the reality can be quite different. When MET is used alongside a
% Cerebus using the Serial port to trigger ephys recording, it is possible
% for human error to corrupt the naming of data, putting the Cerebus and
% MET records our of sync. This tool exists to help repair that damage.
% 
% makrepair( 'nspcheck' , sess )
% 
%   It is possible that a dropped packet or a disruption to File Storage
%   resulted in one or more trials having lost their Cerebus data files.
%   The other consequence is that after the anomoly, Cerebus file names can
%   become out of sync with MET trial ID's. nspcheck is used to determine
%   whether there is a set of Cerebus files for each valid trial ID in the
%   named MET session. sess names a MET session directory.
%   
%   e.g. in a date directory we have MET session directory
%   
%   M123.321.3.Third.Session
%   
%   Assuming that the date directory is the present working directory, as
%   named by pwd, then the command ...
%   
%   makrepair( 'nspcheck' , 'M123.321.3.Third.Session' )
%   
%   Will check that there is a directory called ./nsp/M123_321_3* that
%   contains sets of files named trial_<tid>.<suf> where there is a group
%   of trials (different file type suffixes <suf>) that exists for every
%   trial id <tid> from 1 to the trial_count in the MET session footer.
% 
% 
% makrepair( 'nspincr' , sess , tid1 , tid2 , inc )
%   
%   We may find that the Cerebus NSP trial ID's become out of sync with the
%   MET session directory. This could result from disruption to the File
%   Storage application in Cerebus > Central, throwing the automated
%   incrementation of trial numbers off from that used by MET. nspincr will
%   find a subset of NSP files with trial IDs from tid1 to tid2 and then
%   rename each file such that the trial ID is incremented by inc. tid2 and
%   inc are optional, these default to the last valid trial ID in the MET
%   session for tid2 and a value of +1 for inc. sess is as before in
%   nspcheck.
%   
%   e.g. in a date directory we have MET session dir
%   
%   M123.321.3.Third.Session
%   
%   And NSP session directory
%   
%   nsp/M123_321_3
%   
%   Sadly, the NSP file names go out of sync with the MET trials from trial
%   ID 366. We can repair this by incrementing the NSP file name trial IDs
%   by +1 using
%   
%   makrepair( 'nspincr' , 'M123.321.3.Third.Session' , 366 )
%   
% 
% makrepair( 'fork' , sessfrom , sessto )
% 
%   It is possible that MET advances to a new session but that the user
%   forgets to update the session name in the Cerebus > Central > File
%   Storage application. In this case, trials are added under the previous
%   session name. To repair this, we must create a new nsp sub-directory
%   for the new session. The excess files from the old session are copied
%   over and renamed with a corrected trial index. sessfrom names the full
%   or relative path to the session directory with excess nsp data. sessto
%   names the full or relative path to the session directory that is
%   missing nsp data due to the naming error.
%   
%   e.g. in a date directory we have MET session directories ...
%   
%   M123.321.2.Second.Session --> 100  trials collected
%   M123.321.3.Third.Session  --> 1000 trials collected
%   
%   But in the nsp directory we only have a sub-directory for the second
%   session ...
% 
%   nsp/M123_321_2 --> Contains 1100 .nev files
%   
%   The user forgot to tell Cerebus > Central > File Recording to rename
%   trials 1001 to 1100 for the third session.
%   
%   Run ...
%   
%   sessfrom = 'M123.321.2.Second.Session' ;
%   sessto   = 'M123.321.3.Third.Session'  ;
%   makrepair( 'fork' , sessfrom , sessto )
%   
%   To create directory ...
%   
%   nsp/M123_321_3
%   
%   ... and then move/rename excess trials from session 2 to session 3 ...
%   
%   nsp/M123_321_2/trial_101.* --> nsp/M123_321_3/trial_1.*
%   nsp/M123_321_2/trial_102.* --> nsp/M123_321_3/trial_2.*
%   nsp/M123_321_2/trial_103.* --> nsp/M123_321_3/trial_3.*
%   .
%   .
%   .
%   nsp/M123_321_2/trial_1100.* --> nsp/M123_321_3/trial_1000.*
% 
% Written by Jackson Smith - January 2020 - DPAG , University of Oxford
% 
  
  % Check function string
  if  ~ iscellstr( { funstr } )
    
    error (  'MAK:makrepair:funstr'  ,  ...
      'makrepair: funstr must be a valid character string'  )
    
  end % check funstr
  
  % Obtain MET constants
  [ MC , MCC ] = makmetcon ;
  C.MC = MC ;  C.MCC = MCC ;
  
  % Other constants
  
    % Cerebus nsp sub-directory
    C.nsp = 'nsp' ;
  
  % Run named function
  switch  funstr
    
    % Fork ephys data from one MET session to another
    case  'fork'  ,  mr_fork( C , varargin{ : } )
      
    % Make sure that there is a Cerebus file group for each valid trial ID
    case  'nspcheck'  ,  mr_nspcheck( C , varargin{ : } )
      
    % Increment Cerebus NSP trial ID's to re-set alignment with MET dir
    case  'nspincr'  ,  mr_nspincr( C , varargin{ : } )
    
    % No such function
    otherwise
      
      error (  'MAK:makrepair:invalidfunction'  ,  ...
        'makrepair: Unrecognised function name ''%s'''  ,  funstr  )
      
  end % run named function
  
end % makrepair


%%% Sub-routines %%%

function  mr_fork( C , sessfrom , sessto )
  
  % MET controller constants
  MC = C.MC ;  MCC = C.MCC ;
  
  % Pack into a cell array
  sess = { sessfrom , sessto } ;

  % Guarantee that inputs are character strings
  if  ~ iscellstr( sess )
    
    error (  'MAK:makrepair:fork_cellstr'  ,  [ 'makrepair:fork: ' , ...
      'sessfrom and sessto must both be character strings' ]  )
    
  % Check that both are valid directories
  elseif  ~ exist( sessfrom , 'dir' )  ||  ~ exist( sessto , 'dir' )
    
    error (  'MAK:makrepair:fork_existdir'  ,  [ 'makrepair:fork: ' , ...
      'sessfrom and sessto must both name existing directories' ]  )
    
  end % check input
  
  % Convert to absolute paths
  sess = abpath( sess ) ;
  
  % Break into parent (probably date) directory and MET session directory
  parent = cf( @fileparts , sess ) ;
  
  % Check that session directories are from the same experiment
  if  ~ strcmp( parent{ : } )
    
    error (  'MAK:makrepair:fork_parentdir'  ,  [ 'makrepair:fork: ' , ...
      'sessfrom and sessto have different parent directories' ]  )
    
  end % check same parent dirs
  
  % Collapse to single string
  parent = parent{ 1 } ;
  
  % NSP directory name
  nspdir = fullfile( parent , C.nsp ) ;
  
  % Look for nsp directory
  if  ~ exist( nspdir , 'dir' )
    
    error (  'MAK:makrepair:fork_nspdir'  ,  [ 'makrepair:fork: ' , ...
      'can''t find nsp directory %s' ]  ,  nspdir  )
    
  end % check nsp dir
  
  % Session descriptor and footer files
  sdnam = cf( @( s ) fullfile( s , MCC.SDFNAM  ) , sess ) ;
  footn = cf( @( s ) fullfile( s , MC.SESS.FTR ) , sess ) ;
  
  % Check that all files exist
  for  F = [ sdnam , footn ]  ,  f = F{ 1 } ;
    
    if  ~ exist( f , 'file' )
      error (  'MAK:makrepair:fork_existfile'  ,  [ 'makrepair:fork: ' ,...
      'can''t find file %s' ]  ,  f  )
    end
    
  end % check files exist
  
  % Load session descriptors and footer structs
  sd = cf( @( s ) lv( s , 'sd' ) , sdnam ) ;
  ft = cf( @( s ) lv( s , 'f'  ) , footn ) ;  ft = [ ft{ : } ] ;
  
  % nsp session sub-directory search strings
  nspsub = cf( ...
    @( sd ) fullfile( nspdir , sprintf( '%s_%d_%d*/trial_*' , ...
      sd.subject_id , sd.experiment_id , sd.session_id ) ) , sd ) ;
	
	% Look for nsp trial data file names
  tname = cf( @dir , nspsub ) ;
  
  % We shoud find something in sessionfrom's directory, and nothing in
  % sessionto's
  if  isempty( tname{ 1 } )  ||  ~ isempty( tname{ 2 } )
    
    error (  'MAK:makrepair:fork_nspfiles'  ,  [ 'makrepair:fork: ' ,...
      'files missing from sessionfrom''s nsp dir or sessionto''s ' , ...
        'nsp dir already has data' ]  )
    
  end % check nsp dirs
  
  % Collapse to single set
  tname = tname{ 1 } ;
  
  % Count excess trial groups in sessfrom nsp directory
  [ n , I ] = trialcount( tname , ft( 1 ).trial_count + 1 , ...
    sum( [ ft.trial_count ] ) ) ;
  
  % Make sure we have enough to fork over to the other session
  if  ft( 2 ).trial_count  ~=  n
    
    error (  'MAK:makrepair:fork_mismatch'  ,  [ 'makrepair:fork: ' ,...
      'sessto requires %d nsp trial groups from sessfrom, ' , ...
        'but %d were found' ]  ,  ft( 2 ).trial_count  ,  n  )
    
  end % check number of trials
  
  % sessto nsp sub-directory name , kick off */trial_*
  nspto = regexprep( nspsub{ 2 } , '\*/trial_\*$' , '' ) ;
  
  % Check that directory does not exist
  if  exist( nspto , 'dir' )
    
    error (  'MAK:makrepair:fork_nspto'  ,  [ 'makrepair:fork: ' , ...
      'nsp sub-dir already exists %s' ]  ,  nspto  )
    
  end % check for sessto nsp dir
  
  % Create directory
  fprintf( 'Creating %s\n' , nspto )
  mkdir( nspto )
  
  % Verify directory creation
  if  ~ exist( nspto , 'dir' )
    
    error (  'MAK:makrepair:fork_nspto_fail'  ,  [ 'makrepair:fork: ' , ...
      'failed to create nsp sub-dir %s' ]  ,  nspto  )
    
  end % check sessto nsp dir
  
  % sessfrom's excess files , tid refers to sessto trial id
  for  tid = 1 : n
    
    % sessfrom files
    for  i = I{ tid }
      
      % Full path of file in sessfrom
      nfrom = fullfile( tname( i ).folder , tname( i ).name ) ;
      
      % File type suffix
      suff = regexp( tname( i ).name , '\.(\w{3})$' , 'tokens' ) ;
      suff = suff{ 1 }{ 1 } ;
      
      % Full path of new file name in sessto
      nto = fullfile( nspto , sprintf( 'trial_%d.%s' , tid , suff ) ) ;
      
      % Report
      fprintf( 'Moving %s\n    to %s\n' , nfrom , nto )
      
      % Check that new file name is not already in use
      if  exist( nto , 'file' )
        
        error (  'MAK:makrepair:fork_nto'  ,  [ 'makrepair:fork: ' , ...
          'name already in use %s' ]  ,  nto  )
        
      end % nto already in use
      
      % Simultaneously move and rename the file
      [ s , m ] = movefile( nfrom , nto ) ;
      
      % Error
      if  ~ s
        error( 'MAK:makrepair:fork_movefile' , [ 'makrepair:fork:' , ...
          'failed to move %s to %s\nmovefile error: %s' ] , nfrom , ...
            nto , m )
      end
      
    end % sessfrom files
    
  end % excess files
  
end % mr_fork


function  mr_nspcheck( C , sess )
  
  % Things we need
  [ sess , ~ , trials , ~ , ft , n , I ] = loadstuff( C , sess ) ;
  
  % MET signal struct, map names to codes
  MSIG = C.MC.SIG' ; MSIG = struct( MSIG{ : } ) ;
  
  % Does the count match the number of trials in the MET session?
  if  ft.trial_count  ~=  n
    
    warning( 'MAK:makrepair:nspcheck_count' , [ 'makrepair:nspcheck:' , ...
      '%d trials reported but Cerebus nsp data only found for %d ' , ...
        'trials in MET session %s' ] , ft.trial_count , n , sess )
    
  end % is there nsp trial data for each valid trial id?
  
  % Now we will step through and check whether MET signal events that were
  % recorded by the Cerebus NSP match the number and order in which they
  % occurred. Step through trial IDs.
  for  tid = 1 : ft.trial_count
    
    % Report
    fprintf( 'Trial ID %d of %d\n' , tid , ft.trial_count )
    
    % String version of trial ID
    strtid = sprintf( '%d' , tid ) ;
    
    % Name of trial directory
    tridir = fullfile( sess , C.MC.SESS.TRIAL , strtid ) ;
    
    % Load MET signal data
    msig = load( fullfile( tridir , [ 'metsigs_' , strtid , '.mat' ] ) ) ;
    
    % Load Cerebus event file
    t = trials( I{ tid } ) ;
    
    % Get the one ending in .nev
    t = t( cellfun( @( s ) strcmp( s( end - 2 : end ) , 'nev' ) , ...
      { t.name } ) ) ;
    
    % Can't find it , skip to next trial if it is missing
    if  isempty( t )
      warning( 'MAK:makrepair:nspcheck_nev' , [ 'makrepair:nspcheck:' , ...
        'can''t find .nev file for trial %d of %s' ] , tid , sess )
      continue
    end % check for .nev file
    
    % Load NSP event data
    nev = openNEV( fullfile( t.folder , t.name ) , 'nosave' ) ;
    nev = nev.Data.SerialDigitalIO ;
    
    % Starting index in MET record and NSP record will be the mstart signal
    metmsi = find( msig.sig == MSIG.mstart , 1 , 'first' ) ;
    nevmsi = find( nev.UnparsedData == MSIG.mstart  &  ...
      nev.InsertionReason == 1 , 1 , 'first' ) ;
    
    % Make sure we found mstart signal records
    if  isempty( metmsi )  ||  isempty( nevmsi )
      warning( 'MAK:makrepair:nspcheck_nev' , [ 'makrepair:nspcheck:' , ...
        'can''t find mstart signal for trial %d of %s' ] , tid , sess )
      continue
    end
    
    % Fine, we found mstart, now go to the next event
    metmsi = metmsi  +  1 ;
    nevmsi = nevmsi  +  1 ;
    
    % Number of NSP digin values
    ndigin = numel( nev.UnparsedData ) ;
    
    % NSP record of MET signals
    while  nevmsi <= ndigin
      
      % This is a cargo value or Serial port event
      if  256 <= nev.UnparsedData( nevmsi )  ||  ...
            1 ~= nev.InsertionReason( ( nevmsi ) )
        
        % Advance index and go to next digital event
        nevmsi = nevmsi  +  1 ;
        continue
        
      end % not MET signal code
      
      % We've exceeded the MET session directory's record of the number of
      % MET signals, but not the NSP's record
      if  numel( msig.sig ) < metmsi
        
        warning( 'MAK:makrepair:nspcheck_numsig', ...
          [ 'makrepair:nspcheck:number of MET signals misaligned ' , ...
            'trial %d of %s' ] , tid , sess )
        break
      
      % Check that MET signal code is the same , if not then skip to next
      % trial
      elseif  msig.sig( metmsi )  ~=  nev.UnparsedData( nevmsi )
        
        warning( 'MAK:makrepair:nspcheck_sig', [ 'makrepair:nspcheck:',...
        'MET signal misalignment trial %d of %s' ] , tid , sess )
        break
        
      % Now make sure that the cargo value is the same, unless the cargo is
      % 256 or more
      elseif  nevmsi < ndigin  &&  msig.crg( metmsi ) < 256  &&  ...
          256 <= nev.UnparsedData( nevmsi + 1 )  &&  ...
            msig.crg( metmsi ) ~= ...
              bitshift( nev.UnparsedData( nevmsi + 1 ) , -8 )
        
        warning( 'MAK:makrepair:nspcheck_crg', [ 'makrepair:nspcheck:',...
        'MET cargo misalignment trial %d of %s' ] , tid , sess )
        break
        
      end % check MET signal and cargo values
      
      % Advance MET signal indices
      metmsi = metmsi  +  1 ;
      nevmsi = nevmsi  +  1 ;
      
    end % NSP MET signals
    
  end % trial IDs
    
end % mr_nspcheck


function  mr_nspincr( C , sess , tid1 , tid2 , inc )
  
  [ ~ , ~ , trials , ~ , ft , ~ , I ] = loadstuff( C , sess ) ;
  
  % Make sure that tid1 is within range of valid values
  if  tid1 < 1 || ft.trial_count < tid1
    
    error( 'MAK:makrepair:nspincr_tid1' , [ 'makrepair:nspincr:' , ...
      'tid1 is outside valid range of 1 to %d' ] , ft.trial_count )
    
  end
  
  % tid2 given?
  if  nargin < 4
    
    % No, set default
    tid2 = ft.trial_count ;
    
  % Yes, check valid range
  elseif  tid2 < 1 || ft.trial_count < tid2
    
    error( 'MAK:makrepair:nspincr_tid2' , [ 'makrepair:nspincr:' , ...
      'tid2 is outside valid range of 1 to %d' ] , ft.trial_count )
    
  end
  
  % inc given? If not then set default
  if  nargin < 5  ,  inc = 1 ;  end
  
  % Order of trial depends on whether we are incrementing or decrementing
  % trial IDs. If incrementing then go backwards; otherwise go forwards.
  % This avoids the risk of trying to rename a file to a file name that is
  % already in use.
  if  0 < inc
    
    TID = tid2 : -1 : tid1 ;
    
  else
    
    TID = tid1 : tid2 ;
    
  end
  
  % Trials
  for  tid = TID
    
    % New trial ID string
    tidstr = sprintf( '%d' , tid + inc ) ;
    
    % File types
    for  i = 1 : numel( I{ tid } )
      
      % Get dir info for this file
      t = trials( I{ tid }( i ) ) ;
      
      % Old file name
      old = fullfile( t.folder , t.name ) ;
      
      % Replace old trial ID with the incremented one
      new = fullfile( t.folder , ...
        [ 'trial_' , tidstr , old( end - 3 : end ) ] ) ;
      
      % Check that new file name is unused
      if  exist( new , 'file' )
        
        error( 'MAK:makrepair:nspincr_new' , [ 'makrepair:nspincr:' , ...
          'can''t rename, file name already in use %s' ] , new )
        
      end
      
      % Rename file
      [ status , msg ] = movefile( old , new ) ;
      
      % Failed to rename
      if  ~ status
        
        error( 'MAK:makrepair:nspincr_mv' , [ 'makrepair:nspincr:' , ...
          'failed to rename %s to %s\nmovefile error: %s' ] , ...
            old , new , msg )
        
      end
      
      % Report
      fprintf( '%s renamed to\n\t%s\n' , old , new )
      
    end % types
    
  end % trials
  
end % mr_nspincr


% Takes cell array of strings f naming files and returns full, absolute
% path name of each in cell array a, where a{ i } is the absolute path of
% file f{ i }.
function  a = abpath( f )
  
  % Obtain MATLAB related information about each file, including full path
  s = cf( @what , f ) ;
  
  % Collapse into struct array
  s = [ s{ : } ] ;
  
  % And return absolute paths in cell array
  a = { s.path } ;
  
end % abpath


% In-house version of cellfun in which 'UniformOutput' false is always the
% final set of arguments i.e. cell arrays always returned
function  varargout = cf( varargin )

  varargout = cell( 1 , nargout ) ;
  [ varargout{ : } ] = cellfun( varargin{ : } , 'UniformOutput' , false ) ;

end % cf


% Load and return variable from mat file
function  v = lv( fname , var )
  
  v = load( fname , var ) ;
  v = v.( var ) ;
  
end % load var


% Count number of trials with trial IDs from tid1 to tid2
function  [ n , I ] = trialcount( tname , tid1 , tid2 )
  
  % File name regular expression, use to extract trial id
  rex = 'trial_0*(\d+)\.\w{3}' ;
  
  % NSP file names
  fname = { tname.name } ;
  
  % Extract trial ID strings
  TID = regexp( fname , rex , 'tokens' ) ;
  
  % Look for failure to match trial ID
  i = cellfun( @isempty , TID ) ;
  
  % Put empty string into all empty nested cell arrays
  TID( i ) = repmat( {{ '' }} , 1 , sum( i ) ) ;
  
  % Collapse into a single cell array of strings. Yes, we concatenate
  % twice. See doc regexp for output argument 'out' pertaining to the
  % 'tokens' keyword.
  TID = [ TID{ : } ] ;
  TID = [ TID{ : } ] ;
  
  % And convert to numeric values
  TID = str2double( TID ) ;
  
  % Initialise counter to zero
  n = 0 ;
  
  % Optional indexing data, returns linear indeces of files with same tid
  i = 0 ;
  I = cell( 1 , tid2 - tid1 + 1 ) ;
  
  % trial id's
  for  tid = tid1 : tid2  ,  i = i + 1 ;
    
    % Find linear index vector of non-empty cells
    I{ i } = find( TID == tid ) ;
    
    % Count group of trials
    n = n  +  ~ isempty( I{ i } ) ;
    
  end % trial id's
  
end % trialcount


% Load things that are needed by nspcheck and nspincr
function  [ sess , nspsess , trials , sd , ft , n , I ] = loadstuff( C ,...
  sess )
  
  % MET controller constants
  MC = C.MC ;  MCC = C.MCC ;
  
  % Check that directory exists
  if  ~ exist( sess , 'dir' )
    
    error( 'MAK:makrepair:nspcheck_sess' , [ 'makrepair:nspcheck:' , ...
      'directory not found %s' ] , sess )
    
  end % check dir exists
  
  % Convert to absolute path
  sess = abpath( { sess } ) ;
  sess = sess{ 1 } ;
  
  % Get parent directory
  parent = fileparts( sess ) ;
  
  % Session descriptor and footer files
  sdname = fullfile( sess , MCC.SDFNAM  ) ;
  ftname = fullfile( sess , MC.SESS.FTR ) ;
  
  % Build Cerebus NSP sub-directory name
  nspdir = fullfile( parent , C.nsp ) ;
  
  % Make sure that these things exist
  for  F = { sdname , ftname , nspdir }  ,  f = F{ 1 } ;
    
    if  ~ (  exist( f , 'file' )  ||  exist( f , 'dir' )  )
      
      error( 'MAK:makrepair:nspcheck_exist' , [ 'makrepair:nspcheck:' , ...
      'can''t find %s' ] , f )
    
    end
  end % check existance
  
  % Load session descriptor and footer
  sd = lv( sdname , 'sd' ) ;
  ft = lv( ftname ,  'f' ) ;
  
  % Build MET session's NSP sub-directory
  nspsess = fullfile( nspdir , sprintf( '%s_%d_%d*' , ...
    sd.subject_id , sd.experiment_id , sd.session_id ) ) ;
  
  % Search for directory
  nsp = dir( nspsess ) ;
  
  % Does it exist?
  if  ~ isscalar( nsp )  ||  isempty( nsp )
    
    error( 'MAK:makrepair:nspcheck_nspsub' , [ 'makrepair:nspcheck:' , ...
      'no nsp session dir matching %s' ] , nspsess )
    
  end % look for nsp sess dir
  
  % Get expanded name
  nspsess = fullfile( nsp.folder , nsp.name ) ;
  
  % Look for Cerebus nsp trials
  trials = dir( fullfile( nspsess , 'trial_*' ) ) ;
  
  % Any found?
  if  isempty( trials )
    
    error( 'MAK:makrepair:nspcheck_trials' , [ 'makrepair:nspcheck:' , ...
      'no trial data found in %s' ] , nspsess )
    
  end % look for nsp trial data
  
  % Count each valid trial ID with Cerebus NSP trial data
  [ n , I ] = trialcount( trials , 1 , ft.trial_count ) ;
  
end % loadstuff

