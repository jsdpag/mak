
function  ks = makkeepspk (  ndir ,  sub ,  eid ,  tol_s ,  tol_c ,  ...
                                                                  grpspk  )
% 
% ks = makkeepspk (  ndir ,  sub ,  eid ,  tol_s ,  tol_c ,  grpspk  )
% 
% MET Analysis Kit. Determines which spikes to keep, discarding those that
% appear to result from mere cross-talk between channels. Data is loaded
% from directory ndir (string) for subject sub (string) and experiment eid
% (scalar numeric). A spike is discarded when another unit with higher
% average waveform RMS occurred within tol_s seconds and the two waveforms
% had a Pearson correlation coefficient of at least tol_c. If unspecified,
% tol_s is 4 / 3e4 i.e. 133 micro-seconds or 4 samples at a 30kHz sample
% rate. Likewise, tol_c is 0.5 if unspecified. grpspk is an optional scalar
% logical that defaults to false ; if true then ks is reformatted into an N
% x C cell array with rows indexing trials and columns indexing spike
% clusters in their original order such that ks{ i , j } is the logical
% index variable saying which spikes from trial i on cluster j are kept.
% 
% Returns a cell array ks of size 1 x N where N is the number of trials.
% Each element of ks contains a 1 x M logical vector where M is the number
% of spikes from all electrodes on that trial. The vector is true ( 1 ) for
% kept spikes and false ( 0 ) for discarded spikes.
% 
% The event, spike, and manual spike files are required for each
% experiment.
% 
% Written by Jackson Smith - July 2019 - DPAG, University of Oxford
% 
 

  %%% Check Input %%%
  
  % Number of inputs and outputs
   narginchk (  3  ,  6  )
  nargoutchk (  0  ,  1  )
  
  % No time tolerance given so use default
  if  nargin  <  4  ||  isempty (  tol_s  )  ,  tol_s = 4  /  3e4 ;  end
  
  % No correlation tolerance given so use default
  if  nargin  <  5  ||  isempty (  tol_c  )  ,  tol_c = 0.5 ;  end
  
  % Default grpspk is false
  if  nargin  <  6  ||  isempty (  grpspk  )  ,  grpspk = false ;  end
  
  % ndir is string
  if  ~ ischar (  ndir  )  ||  ~ isrow (  ndir  )
    
    error (  'MAK:makkeepspk:ndirstr'  ,  ...
      'makkeepspk: ndir must be a string'  )
    
  % Does directory exist?
  elseif  ~ exist (  ndir  ,  'dir'  )
    
    error (  'MAK:makkeepspk:nodir'  ,  ...
      'makkeepspk: can''t find directory %s'  ,  ndir  )
    
  % sub must be a string
  elseif  ~ ischar (  sub  )  ||  ~ isrow (  sub  )
    
    error (  'MAK:makkeepspk:substr'  ,  ...
      'makkeepspk: sub must be a string'  )
    
  % eid must be a positive non-zero integer
  elseif  ~ isscalar (  eid  )  ||  ~ isreal (  eid  )  ||  ...
      ~ isfinite (  eid  )  ||  mod (  eid  ,  1  )  ||  eid  <=  0
    
    error (  'MAK:makkeepspk:eid'  ,  ...
      'makkeepspk: eid must be a scalar, finite integer of 1 or more'  )
    
  % tol_s must be scalar, real, finite and non-negative
  elseif  ~ isscalar (  tol_s  )  ||  ~ isreal (  tol_s  )  ||  ...
      ~ isfinite (  tol_s  )  ||  tol_s  <  0
    
    error (  'MAK:makkeepspk:tol_s'  ,  [ 'makkeepspk: tol_s must be ' ,...
      'a scalar, real, finite and non-negative' ]  )
    
  % tol_c must be scalar, real, finite and in the range [ -1 , +1 ]
  elseif  ~ isscalar (  tol_c  )  ||  ~ isreal (  tol_c  )  ||  ...
      ~ isfinite (  tol_c  )  ||  tol_c  <  -1  ||  tol_c  >  +1
    
    error (  'MAK:makkeepspk:tol_c'  ,  [ 'makkeepspk: tol_c must be ' ,...
      'a scalar, real, finite and between -1 and +1, inclusive' ]  )
    
  % grpspk must be a scalar logical
  elseif  ~ isscalar (  grpspk  )  ||  ~ islogical (  grpspk  )
    
    error (  'MAK:makkeepspk:eid'  ,  ...
      'makkeepspk: grpspk must be a scalar logical'  )
    
  end % check input
  
  % File base name
  fbase = sprintf (  '%s.%d.*'  ,  sub  ,  eid  ) ;
  
  % Report
  fprintf (  'Loading %s\n'  ,  fbase  )
  
  % Look for event, spike, and manual spike files
  df = dir ( fullfile(  ndir  ,  [ fbase , 'mat' ]  ) ) ;
  
  % There are at least three files expected
  if  numel (  df  )  <  3
    
    error (  'MAK:makkeepspk:numfiles'  ,  ...
      'makkeepspk: 3 pre-processing files expected but %d found'  ,  ...
        numel( df )  )
    
  end % 3 files expected
  
  % Get spike automated and manual spike sorting file index
  af = ~ cellfun (  @isempty ,  regexp( { df.name } , '.spksort.mat$' )  );
  mf = ~ cellfun (  @isempty ,  regexp( { df.name } ,  '.manual.mat$' )  );
  
  % One of each
  if  sum (  af  )  ~=  1  ||  sum (  mf  )  ~=  1
    
    error (  'MAK:makkeepspk:numspksortfiles'  ,  [ 'makkeepspk: ' ,...
      'one .spksort.mat and one .manual.mat files expected' ]  )
    
  end % one of each
  
  % Get file names
  af = df( af ).name ;
  mf = df( mf ).name ;
  df = strrep (  af ,  '.spksort' ,  ''  ) ;
  
  % Load data
  d = load (  fullfile( ndir , df )  ,  'd'  ) ;  d = d.d ;
  s = load (  fullfile( ndir , af )  ,  's'  ) ;  s = s.s ;
  m = load (  fullfile( ndir , mf )  ,  'm'  ) ;  m = m.m ;
  
  % Do quick check on loaded data
  if  ~ isstruct (  d  )  ||  ~ isstruct (  s  )  ||  ~ isstruct (  m  )
    
    error (  'MAK:makkeepspk:nostructs'  ,  [ 'makkeepspk: ' ,...
      '%s.mat files must contain structs d, s, and m' ]  ,  fbase  )
    
  % Scalar structs expected
  elseif  ~ isscalar (  d  )  ||  ~ isscalar (  s  )  ||  ...
      ~ isscalar (  m  )
    
    error (  'MAK:makkeepspk:scalarstructs'  ,  [ 'makkeepspk: ' ,...
      'structs d, s, and m must all be scalar' ]  )
    
  % Event data must contain certain key fields
  elseif  ~ all ( isfield(  d  ,  ...
      { 'preppar' , 'subject_id' , 'start' , 'spike' }  ) )
    
    error (  'MAK:makkeepspk:dfields'  ,  [ 'makkeepspk: ' ,...
      'd appears to have incorrect fields' ]  )
    
  % Automated spike sorting fields
  elseif  ~ all ( isfield(  s  ,  ...
      { 'wave_raw' , 'wave_aligned' , 'pca' , 'E_init' }  ) )
    
    error (  'MAK:makkeepspk:sfields'  ,  [ 'makkeepspk: ' ,...
      's appears to have incorrect fields' ]  )
    
	% Manual spike sorting fields
  elseif  ~ all ( isfield(  m  ,  ...
      { 'cutoff' , 'clustmap' , 'wave_avg' , 'clustind' }  ) )
    
    error (  'MAK:makkeepspk:mfields'  ,  [ 'makkeepspk: ' ,...
      'm appears to have incorrect fields' ]  )
    
  end % check structs
  
  
  %%% Preparation %%%
  
  % cellfun UniformOutput false name/value pair
  UF = {  'UniformOutput'  ,  false  } ;
  
  % Cast spike times to double
  spike = cellfun (  @double  ,  d.spike.time  ,  UF{ : }  ) ;
  
  % Electrod ID on each spike
  EID = d.spike.electrode ;
  
  % Linear index to electrode id map
  eidmap = d.electrodes ;
  
  % Cluster ID on each spike
  CID = m.clustind ;
  
  % Grab certain key data from spike sort and manual cluster files
  % including raw waveforms cast to single floating point, ...
  W = cellfun (  @single  ,  s.wave_raw  ,  UF{ : }  ) ;
  
  % ... RMS squared vector, ...
  rms = [ m.rms{ : } ] ;
  
  % ... Electrode ID vector in register with rms, ...
  eid = arrayfun (  @( c , t ) repmat( t , 1 , c )  ,  m.numclust'  ,  ...
    eidmap  ,  UF{ : }  ) ;
  eid = [ eid{ : } ] ;
  
  % ... Cluster ID vector in register with rms
  cid = arrayfun (  @( c ) 1 : c  ,  m.numclust  ,  UF{ : }  ) ;
  cid = [ cid{ : } ] ;
  
  % Sort rms from largest to smallest and order eid and cid accordingly
  [ rms , i ] = sort (  rms  ,  'descend'  ) ;
  eid = eid( i ) ;
  cid = cid( i ) ;
  
    % We will do spike grouping , so make sure that we can restore the
    % original order of spike clusters
    if  grpspk  ,  [ ~ , restore ] = sort (  i  ) ;  end
  
  % Number of clusters
  n.clust = numel (  rms  ) ;
  
  % Number of samples per waveform
  n.samp = unique ( cellfun(  @( w ) size( w , 1 )  ,  W  ) ) ;
  
    if  ~ isscalar (  n.samp  )
      error (  'MAK:makkeepspk:numsamp'  ,  ...
        'makkeepspk: Mixed number of samples per spike waveform'  )
    end
    
	% Number of trials
  n.numtrials = d.numtrials ;
  
  % Number of electrodes
  n.numtrodes = d.numtrodes ;
  
  % Waveform index. W is a cell array with all waveforms from electrode e
  % in W{ e }. wi( e ) is the index of the spike waveform in W{ e } from
  % the previous trial. Add 1 and we get the first spike from the next
  % trial, hence we initialise to zero.
  wi = zeros ( size(  W  ) ) ;
  
  % Allocate output cell array
  ks = cell (  1  ,  d.numtrials  ) ;
  
  % Clear memory of unused loaded data
  clear  d s m
  
  
  %%% Find kept spikes %%%
  
  % Report
  fprintf (  'Finding kept spikes\n'  )
  
  % Trials
  for  t = 1 : n.numtrials
    
    % Report every tenth trial
    if  mod (  t ,  10  )  ==  0
      
      fprintf (  'Trial %d of %d\n'  ,  t  ,  n.numtrials  )
      
    end % report
    
    % Total number of spikes in this trial
    n.totspk = numel (  spike{ t }  ) ;
    
    % Allocate logical vector marking kept spikes , initialise true ( 1 )
    ks{ t } = true (  1  ,  n.totspk  ) ;
    
    % Temporal tolerance bands around each spike
    tmin = spike{ t }  -  tol_s ;
    tmax = spike{ t }  +  tol_s ;
    
    % Collect waveforms from this trial
    w = zeros (  n.samp  ,  n.totspk  ,  'single'  ) ;
    
    % Electrodes
    for  e = 1 : n.numtrodes
      
      % Spikes from this electrode
      i = EID{ t }  ==  eidmap( e ) ;
      
      % Number of spikes
      n.spk = sum (  i  ) ;
      
      % Copy waveform
      w( : , i ) = W{ e }( : , wi( e ) + ( 1 : n.spk ) ) ;
      
      % Advance waveform index on this electrode
      wi( e ) = wi( e )  +  n.spk ;
      
    end % trodes
    
    % Spike clusters from biggest RMS to smallest
    for  c = 1 : n.clust
      
      % Find kept spikes from this cluster in this trial
      i = ks{ t }  &  EID{ t } == eid( c )  &  CID{ t } == cid( c ) ;
      
      % Linear index
      i = find (  i  ) ;
      
      % Spikes from this cluster
      n.spk = numel (  i  ) ;
      
      % Spikes
      for  s = 1 : n.spk
        
        % Find spikes within temporal tolerance , synchronous spikes
        j = tmin <= spike{ t }( i( s ) )  &  ...
            tmax >= spike{ t }( i( s ) ) ;
          
        % Exclude the spike in question
        j( i( s ) ) = 0 ;
        
        % Calculate Pearson correlation coefficient of synchronous spikes
        rho = corr (  w( : , i( s ) )  ,  w( : , j )  ) ;
        
        % Discard any synchronous spikes with waveforms that are too
        % similar to the spike in question
        ks{ t }( j ) = tol_c  >  rho ;
        
      end % spikes
      
    end % spike clusters
    
  end % trials
  
  
  %%% Group spikes %%%
  
  % No spike grouping , return ks as it is
  if  ~ grpspk  ,  return  ,  end
  
  % Restore original order of electrod ID's and cluster ID's
  eid = eid( restore ) ;
  cid = cid( restore ) ;
  
  % Spike grouping cell array
  g = cell (  n.numtrials  ,  n.clust  ) ;
  
  % Trials
  for  t = 1 : n.numtrials
    
    % Spike clusters
    for  c = 1 : n.clust
      
      % Locate spikes on this trial from this cluster
      i = EID{ t } == eid( c )  &  CID{ t } == cid( c ) ;
      
      % Pull out logical vector identifying kept spikes for this cluster on
      % this trial
      g{ t , c } = ks{ t }( i ) ;
      
    end % spk clust
    
  end % trials
  
  % Replace output argument
  ks = g ;
  
  
end % makkeepspk

