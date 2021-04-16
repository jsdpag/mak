
function  [ f , r , p , e , c , n , cl , a ] = ...
  makcrosstalk (  ndir ,  sub ,  eid ,  tol ,  clusters  )
% 
% [ f , r , p , e , c , n ] = makcrosstalk (  ndir ,  sub ,  eid ,  tol  )
% [ ... , cl , a ] = makcrosstalk( ... , '-clusters' )
% 
% MET Analysis Kit. Computes measures of cross-talk between separate
% channels that were recorded simultaneously from subject sub (string)
% in experiment eid (integer, positive, non-zero). Pre-processing event,
% automated and manual spike-sorting files for this subject and experiment
% are found in the directory named ndir.
% 
% Large numbers of synchronous spikes with similar waveforms may indicate
% cross-talk between channels. In a linear probe, this will be expected for
% nearby channels. But in a Utah array with large (e.g. 400 micrometer)
% spacing, it is a sign of equipment or software failure and must be
% accounted for. Hence, spikes from separate electrodes with identical
% spike times are sought ; only spikes that weren't discarded during
% sorting are considered. If spikes a and b with times t( a ) and t( b )
% are within tol seconds of each other i.e. abs( t( a ) - t( b ) ) <= tol
% then they are considered synchronous.  If tol is not provided (3 input
% args or tol is empty i.e. [ ]) then a default value of 4 / 3e4 seconds is
% used i.e. 4 samples at 30,000Hz sampling rate.
% 
% The fraction of spikes with identical spike times (within tolerance range
% +/-tol) is returned in square matrix f( : , : , 1 ), where f( i , j , 1 )
% is the fraction of spikes from the ith electrode that have identical
% spike times on the jth electrode. Matrix f( : , : , 2 ) contains the
% total number of spikes from the ith electrode. These values can be used
% to deduce the underlying binomial distribution.
% 
% The median Pearson correlation coefficient of spike waveforms with
% identical spike times is returned in square matrix r( : , : , 1 ), where
% r( i , j , 1 ) is the median correlation of synchronous spikes between
% the ith and jth electrodes i and j. Matrices r( : , : , 2 : 3 ) contain
% the 2.5 and 97.5 percentiles of the distribution of correlations for
% synchronous spikes.
% 
% Lastly, a set of two-sample tests are used to compare the distribution of
% correlations from synchronous spikes against the correlations of all
% non-synchronous spikes and the resulting p-values are returned in square
% matrix p. p( i , j , : ) are the p-values for tests comparing the
% correlations of synchronous spikes from the ith and jth electrodes. A
% p-value of one is returned in any case when there were no synchronouse
% spikes, or there was just one spike on either electrode. Tests are
% indexed along the 3rd dimension i.e. the layer or page. T-test p-values
% are in layer 1, Wilcoxon rank-sum p-values are in layer 2, and Kolmogorov
% -Smirnov p-values are in layer 3.
% 
% e and c contain the electrode ID's and corresponding channel numbers for
% each row/column of f, r, and p. n is the number of kept electrodes.
% 
% If no output arguments are requested then the output variables will be
% written to a file with format < sub >.< eid >[.< tags >].crosstalk.mat,
% where tags are the same as for the makprep output files. Will attempt to
% remove write permissions from all users.
% 
% It may be desirable instead to compute crosstalk measures between all
% spike clusters. This can be done by adding the optional flag '-clusters'
% as the last input argument. In this case, output arg n is the number of
% kept spike clusters according to the manual spike sorting data. Rows and
% columns of f, r, and p are indexed across spike clusters rather than
% electrodes. e and c still name the electrode ID and channel ID of each
% row/column of f, r, and p but there can now be repeat values for
% electrodes with multiple spike clusters. Output arguments cl and a are
% also returned. cl is the cluster ID, in register with e and c. a is the
% amplitude of the average spike waveform for each spike cluster (max-min),
% also in register with e and c.
% 
% Uses Matlab's Parallel Computing Toolbox.
% 
% At present, function assumes that only one set of pre-processing files
% exists for experiment eid. This may not be the case if different types of
% session were run under experiment eid. For now, a workaround is to create
% symbolic links of the desired files in a separate working directory.
% 
% Written by Jackson Smith - July 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check Input %%%
  
  % Number of inputs and outputs
   narginchk (  3  ,  5  )
  nargoutchk (  0  ,  3  )
  
  % No tolerance given so use default
  if  nargin  <  4  ||  isempty( tol )  ,  tol = 4  /  3e4 ;  end
  
  % ndir is string
  if  ~ ischar (  ndir  )  ||  ~ isrow (  ndir  )
    
    error (  'MAK:makcrosstalk:ndirstr'  ,  ...
      'makcrosstalk: ndir must be a string'  )
    
  % Does directory exist?
  elseif  ~ exist (  ndir  ,  'dir'  )
    
    error (  'MAK:makcrosstalk:nodir'  ,  ...
      'makcrosstalk: can''t find directory %s'  ,  ndir  )
    
  % sub must be a string
  elseif  ~ ischar (  sub  )  ||  ~ isrow (  sub  )
    
    error (  'MAK:makcrosstalk:substr'  ,  ...
      'makcrosstalk: sub must be a string'  )
    
  % eid must be a positive non-zero integer
  elseif  ~ isscalar (  eid  )  ||  ~ isreal (  eid  )  ||  ...
      ~ isfinite (  eid  )  ||  mod (  eid  ,  1  )  ||  eid  <=  0
    
    error (  'MAK:makcrosstalk:eid'  ,  ...
      'makcrosstalk: eid must be a scalar, finite integer of 1 or more'  )
    
  % tol must be scalar, real, finite and non-negative
  elseif  ~ isscalar (  tol  )  ||  ~ isreal (  tol  )  ||  ...
      ~ isfinite (  tol  )  ||  tol  <  0
    
    error (  'MAK:makcrosstalk:tol'  ,  [ 'makcrosstalk: tol must be ' ,...
      'a scalar, real, finite and non-negative' ]  )
    
  end % check input
  
  % -clusters flag?
  if  nargin  <  5
    
    % Nope, lower internal cluster flag
    cluflg = false ;
    
  % Yep! Raise internal cluster flag
  elseif  strcmp( clusters , '-clusters' )
    
    cluflg = true ;
    
  else
    
    error( 'MAK:makcrosstalk:clusters' , ...
      'makcrosstalk:invalid input arg 5, expected it to be ''-clusters''' )
    
  end
  
  % File base name
  fbase = sprintf (  '%s.%d.*'  ,  sub  ,  eid  ) ;
  
  % Report
  fprintf (  'Loading %s\n'  ,  fbase  )
  
  % Look for event, spike, and manual spike files
  df = dir ( fullfile(  ndir  ,  [ fbase , 'mat' ]  ) ) ;
  
  % There are three files expected
  if  any ( ~ cellfun(  @isempty ,  ...
      regexp( { df.name } , '.crosstalk.mat$' )  ) )
    
    error (  'MAK:makcrosstalk:crosstalkfound'  ,  [ 'makcrosstalk: ' ,...
      'a %s.crosstalk.mat file found' ]  ,  fbase  )
    
  elseif  numel (  df  )  ~=  3
    
    error (  'MAK:makcrosstalk:numfiles'  ,  ...
      'makcrosstalk: 3 pre-processing files expected but %d found'  ,  ...
        numel( df )  )
    
  end % 3 files expected
  
  % Get spike automated and manual spike sorting file index
  af = ~ cellfun (  @isempty ,  regexp( { df.name } , '.spksort.mat$' )  );
  mf = ~ cellfun (  @isempty ,  regexp( { df.name } ,  '.manual.mat$' )  );
  
  % One of each
  if  sum (  af  )  ~=  1  ||  sum (  mf  )  ~=  1
    
    error (  'MAK:makcrosstalk:numspksortfiles'  ,  [ 'makcrosstalk: ' ,...
      'one .spksort.mat and one .manual.mat files expected' ]  )
    
  end % one of each
  
  % Get file names
  i = ~( af | mf ) ;
  af = df( af ).name ;
  mf = df( mf ).name ;
  df = df( i ).name ;
  
  % Output file name
  fname = regexprep (  df  ,  '.mat$'  ,  '.crosstalk.mat'  ) ;
  fname = fullfile (  ndir  ,  fname  ) ;
  
  % Output file already exists
  if  exist (  fname  ,  'file'  )
    
    error (  'MAK:makcrosstalk:fileexists'  ,  [ 'makcrosstalk: ' ,...
      'output file already exists: %s' ]  ,  fname  )
    
  end
  
  % Load data
  d = load (  fullfile( ndir , df )  ,  'd'  ) ;  d = d.d ;
  s = load (  fullfile( ndir , af )  ,  's'  ) ;  s = s.s ;
  m = load (  fullfile( ndir , mf )  ,  'm'  ) ;  m = m.m ;
  
  % Do quick check on loaded data
  if  ~ isstruct (  d  )  ||  ~ isstruct (  s  )  ||  ~ isstruct (  m  )
    
    error (  'MAK:makcrosstalk:nostructs'  ,  [ 'makcrosstalk: ' ,...
      '%s.mat files must contain structs d, s, and m' ]  ,  fbase  )
    
  % Scalar structs expected
  elseif  ~ isscalar (  d  )  ||  ~ isscalar (  s  )  ||  ...
      ~ isscalar (  m  )
    
    error (  'MAK:makcrosstalk:scalarstructs'  ,  [ 'makcrosstalk: ' ,...
      'structs d, s, and m must all be scalar' ]  )
    
  % Event data must contain certain key fields
  elseif  ~ all ( isfield(  d  ,  ...
      { 'preppar' , 'subject_id' , 'start' , 'spike' }  ) )
    
    error (  'MAK:makcrosstalk:dfields'  ,  [ 'makcrosstalk: ' ,...
      'd appears to have incorrect fields' ]  )
    
  % Automated spike sorting fields
  elseif  ~ all ( isfield(  s  ,  ...
      { 'wave_raw' , 'wave_aligned' , 'pca' , 'E_init' }  ) )
    
    error (  'MAK:makcrosstalk:sfields'  ,  [ 'makcrosstalk: ' ,...
      's appears to have incorrect fields' ]  )
    
	% Manual spike sorting fields
  elseif  ~ all ( isfield(  m  ,  ...
      { 'cutoff' , 'clustmap' , 'wave_avg' , 'clustind' }  ) )
    
    error (  'MAK:makcrosstalk:mfields'  ,  [ 'makcrosstalk: ' ,...
      'm appears to have incorrect fields' ]  )
    
  end % check structs
  
  % cellfun UniformOutput false name/value pair
  UF = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Preparation %%%
  
  % Report
  fprintf (  '  Preparation\n'  )
  
  % Get kept electrode id's
  e = d.electrodes ;
  
  % Calculate crosstalk between spike clusters
  if  cluflg
    
    % Find ordered list of cluster ID's for each electrode
    cl = arrayfun( @( n ) 1 : n , m.numclust' , UF{ : } ) ;
    cl = [ cl{ : } ] ;
    
    % Repeat electrode ID for each corresponding spike cluster
    e = arrayfun( @( e , n ) repmat( e , 1 , n ) , e , m.numclust' , ...
      UF{ : } ) ;
    e = [ e{ : } ] ;
    
    % Cluster Id comparison function
    fcl = @eq ;
    
  % Crosstalk between electrodes
  else
    
    % We need a cluster vector for the data collection for-loop just below
    cl = zeros( size( e ) ) ;
    
    % Cluster id comparison function
    fcl = @lt ;
    
  end % spk clust
  
  % Get channel numbers for each electrode ID
  c = d.elec2chan(  e  ) ;
  
  % Number of electrodes
  Nt = numel( e ) ;
  n = Nt ;
  
  % Number of unique electrode/cluster pairs , including auto-correlations
  Np = (  Nt ^ 2  -  Nt  )  /  2  +  Nt ;
  
  % Convert trial start times from numeric to cell vector
  d.start = num2cell (  d.start  ) ;
  
  % Get lists of spike times and kept waveforms for each kept electrode
  time = cell (  Nt  ,  1  ) ;
  wave = cell (  Nt  ,  1  ) ;
  
  % Electrodes/spike clusters
  for  i = 1 : Nt
    
    % Electrode id
    eid =  e( i ) ;
    
      % Electrode index
      eind = d.electrodes == eid ;
    
    % Cluster id
    cid = cl( i ) ;
    
    % Kept spikes
    kept = cellfun(  @( c , e ) fcl( cid , c( eid == e ) )  , ...
      m.clustind  ,  d.spike.electrode  ,  UF{ : }  ) ;
    
    % Kept spike times , in PTB i.e. absolute time
    time{ i } = cell2mat(  cellfun( ...
      @( s, t, e, c ) s  +  double( t( eid == e & fcl( cid , c ) ) ) , ...
        d.start,  d.spike.time,  d.spike.electrode,  m.clustind,  UF{:} ));
      
    % Kept waveforms
    wave{ i } = s.wave_raw{ eind }(  :  ,  [ kept{ : } ]  ) ;
    
  end % electrodes
  
  % Find amplitudes of spike clusters
  a = cellfun( @( w ) mean( double( w ) , 2 ) , wave , UF{ : } ) ;
  a = cellfun( @( a ) max( a ) - min( a ) , a )' ;
  
  % Get rid of d, s, and m to clear memory. They are no longer needed.
  clear  d  s  m
  
  
  %%% Compute measures of cross-talk %%%
  
  % Report
  fprintf (  '  Computation\n'  )
  
  % Allocate output memory
  f = zeros (  Nt  ,  Nt  ,  2  ) ;
  r = zeros (  Nt  ,  Nt  ,  3  ) ;
  p =  ones (  Nt  ,  Nt  ,  3  ) ;
  
  % Number of spikes
  Ns = [ 0 , 0 ] ;
  
  % Pair counter
  k = 0 ;
  
  % Electrode pairs
  for  i = 1 : Nt
    
    % ith electrode kept spike times and waveforms
    ti = time{ i } ;
    wi = double (  wave{ i }  ) ;
    Ns( 1 ) = numel (  ti  ) ;
    
    for  j = i : Nt
    
      % Increment unique pairs
      k( 1 ) = k  +  1 ;

      % Report
      fprintf (  '  [ %d , %d ] %d of %d ( %0.3f%% )\n'  ,  ...
        i  ,  j  ,  k  ,  Np  ,  100 * k  /  Np  )

      % jth electrode kept spike times and waveforms
      tj = time{ j } ;
      wj = double (  wave{ j }  ) ;
      Ns( 2 ) = numel (  tj  ) ;

      % Find synchronous spikes on both electrodes
      [ si , sj ] = ismembertol (  ti ,  tj ,  tol ,  'DataScale' ,  1  ) ;

        % Keep indices in j for synchronous spikes
        sj = sj( si ) ;

        % Convert i to linear indices
        si = find (  si  ) ;

      % Spikes from ith electrode
      f( i , j , 2 ) = Ns( 1 ) ;
      
      % Spikes from the jth electrode
      f( j , i , 2 ) = Ns( 2 ) ;

      % Synchronouse spikes detected
      if  ~ isempty (  si  )

        % Fraction of spikes on ith electrode that are synchronous with
        % those from jth electrode
        f( i , j , 1 ) = numel (  si  )  /  Ns( 1 ) ;
        
        % Same again from jth electrode's perspective
        f( j , i , 1 ) = numel (  sj  )  /  Ns( 2 ) ;

        % Syncronous spike waveforms
        swi = wi( : , si ) ;
        swj = wj( : , sj ) ;
        
        % Guarantee fresh sliced output variables , otherwise we carry over
        % junk from the last electrode pair
        clear  pcsync  pcnsyn

        % Compute Pearson correlation of synchronous spikes
        parfor  s = 1 : numel (  si  )

          pcsync( s ) = corr (  swi( : , s )  ,  swj( : , s )  ) ;

        end % corr sync spks

        % Similarity of waveforms. Median correlation with 2.5 and 97.5
        % percentiles.
        r( i , j , : ) = prctile (  pcsync  ,  [ 50 , 2.5 , 97.5 ]  ) ;
        
        % Copy over to other side of diagonal
        r( j , i , : ) = r( i , j , : ) ;

        % Skip to next electrode pair if there is just one spike on either
        % electrode
        if  any (  Ns  ==  1  )  ,  continue  ,  end

        % Sample non-synchronous waveform correlations , with replacement
        parfor  s = 1 : max (  1e3  ,  numel( si )  )

          % Sample any two non-synchronous spikes
          x = ceil( rand( 1 , 2 ) .* Ns ) ;

          while  find( x( 1 ) == si )  ==  find( x( 2 ) == sj )
            x( : ) = ceil (  rand( 1 , 2 )  .*  Ns  ) ;
          end

          % Waveform correlation
          pcnsyn( s ) = corr (  wi( : , x( 1 ) )  ,  ...
                                wj( : , x( 2 ) )  ) ; %#ok

        end % corr non-sync spk sample

        % See whether correlation value of waveforms is stronger for
        % synchronous vs non-synchronous spikes
        [ ~ , p( i , j , 1 ) ] =  ttest2 (  pcsync  ,  pcnsyn  ) ;
              p( i , j , 2 )   = ranksum (  pcsync  ,  pcnsyn  ) ;
        [ ~ , p( i , j , 3 ) ] = kstest2 (  pcsync  ,  pcnsyn  ) ;
        
        % Copy to other side of diagonal
        p( j , i , : ) = p( i , j , : ) ;

      end % fraction of synchronous spikes
    
    end % jth electrode
  end % ith electrode
  
  % Output arguments requested , quit now
  if  nargout  ,  return  ,  end
  
  % Save data
  save (  fname ,  'f' ,  'r' ,  'p' ,  'e' ,  'c' ,  'n' ,  'cl' ,  'a'  )
  
  % Command string
  cstr = sprintf (  'chmod a-w %s'  ,  fname  ) ;
  
  % Remove all write privileges
  [  stat  ,  cmd  ] = system (  cstr  ) ;
  
  % Error
  if  stat
    
    error (  'MAK:makcrosstalk:chmod'  ,  [ 'makcrosstalk: failed to ' ,...
      'remove write permissions from %s: %s' ]  ,  fname  ,  cmd  )
    
  end
  
  
end % makcrosstalk

