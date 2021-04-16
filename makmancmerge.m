
function  varargout = makmancmerge (  varargin  )
% 
% m = makmancmerge ( d , s ). Returns the results of manual cluster
%   merging in m based on event data d and automated spike sorting s, both
%   of which are structs returned from makprep.
% 
% m = makmancmerge ( prep ). Reads d and s from files with base name prep,
%   withough .mat suffix. [ prep , '.mat' ] will contain d and [ prep ,
%   '.spksort.mat' ] will contain s.
% 
% makmancmerge ( prep ). Instead of returning m, it is written to a file
%   called [ prep , '.manual.mat' ]. Write permissions are automatically
%   removed.
% 
% makmancmerge ( prep , d , s ). d and s are used for cluster merging, but
%   the result is automatically saved to file, as above.
% 
% MET Analysis Kit, pre-processing. Steps through all electrodes in a data
% set, allowing manual cluster merging to be done for each one. The results
% are gathered together and a final cluster index is assigned to each
% spike.
% 
% The final cluster index is zero for rejected spikes. It is then 1, 2, 3,
% and so on for spikes with increasing root-mean-squared. Hence, spikes in
% cluster 1 will have a lower RMS than those in cluster 2.
%   
% 
% Output
% 
% m - Struct containing the results of manual cluster merging, has fields:
%   
%   .cutoff - Double floating point vector of connection-strength cutoff
%     points for each electrode. These are used to determine when to stop
%     automatic spike cluster merging.
%   
%   .mergers - cell column vector of manual cluster mergers for each
%     electrode. Each element contains a 2 x N uint8 matrix listing N
%     manual cluster mergers. The ith column lists the ith merger between
%     two clusters. Rows 1 contains the low-numbered cluster and row 2
%     contains the high-numbered cluster of the merger. Columns are ordered
%     chronologically so that merger i occurred before i + 1. Initial
%     cluster index numbers are used.
% 
%   .E - cell column vector of interface energy matrices for each
%     electrode. Each element has an N x N double matrix of the final
%     interface energy matrix following all cluster mergers starting from N
%     initial clusters (see makprep). There is a row and column for each
%     initial cluster.
%
%   .clustmap - cell column vector of mappings from initial to final
%     cluster indices for each electrode. Each element has a 1 x N uint8
%     vector mapping the initial cluster index to the final index value.
%     For the ith electrode and jth initial cluster, .clustmap{ i }( j ) is
%     the final cluster value used in .clustind. The initial cluster
%     indices are those used in s (see makprep). Zero-value place holders
%     are put for clusters that were merged into other clusters, and those
%     that were manually rejected. Use these index values to refer to spike
%     clusters in all following fields.
%   
%   .numclust - uint8 column vector. The number of clusters that remain
%     after automated and manual merging on each electrode, not including
%     rejected spikes.
%   
%   .wave_avg - cell column vector of average spike cluster waveforms for
%     each electrode. Each element has an N x M single floating point
%     matrix of M average waveforms each with N samples. The average
%     waveform on the ith electrode for the jth spike cluster is
%     .wave_avg{ i }( : , j ).
%   
%   .wave_var - cell column vector of spike cluster waveform variance for
%     each electrode. Each element has an N x M single floating point
%     matrix of waveform variance for M clusters, each with N samples per
%     waveform. The waveform variance of the ith electrode and jth spike
%     cluster is .wave_var{ i }( : , j ).
%   
%   .rms - cell column vector of root-mean-squared for each cluster on each
%     electrode. Each element contains a 1 x N single floating point vector
%     of the RMS for each of N spike clusters on the electrode. The RMS for
%     the ith electrode and jth spike cluster is .rms{ i }( j ).
%   
%   .count - cell column vector of spike counts per cluster and trial for
%     each electrode. Each element has an N x M uint32 matrix of the spike
%     count for each of N clusters on each of M trials.
% 
%   .clustind - 1 x N cell vector contains the final cluster index for all
%     spikes on each of N trials. Each element contains a uint8 vector of
%     the final cluster index for each spike. For a given trial , i ,
%     .clustind{ i } groups spikes from all electrodes. The data is in
%     register with d.spike.electrode and d.spike.time.
% 
% 
% To analyse the data, one should use d and m together. The majority of
% data will be stored in d, while spike cluster identities are taken from
% m. For example, to get the time stamps of all spikes from trial i, on
% electrode j, in cluster k we would execute:
%   
%   s = d.spike.electrode{ i } == j  &  m.clustind{ i } == k ;
%   t = d.spike.time{ i }( s ) ;
% 
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input arguments %%%
  
  % Minimum/maximum number of input and output arguments
  nargoutchk (  0  ,  1  )
  narginchk  (  1  ,  3  )
  
  % Index prior to first input struct
  i = 0 ;
  
  % Default file names
  fdata = '' ;
  fsort = '' ;
   fout = '' ;
   
  % Default data structs
  d = [] ;
  s = [] ;
  
  % One or three inputs
  if  any(  nargin  ==  [ 1 , 3 ]  )
    
    % First arg must be prep
    if  ~ (  isvector( varargin{ 1 } )  &&  ischar( varargin{ 1 } )  )
    
      error (  'MAK:makmancmerge:prep'  ,  [ 'makmancmerge: input arg ' , ...
        'prep must be a string' ]  )
      
    % Look for .mat suffix
    elseif  ~ isempty ( regexp(  varargin{ 1 }  ,  '.mat$'  ,  'once'  ) )
      
      error (  'MAK:makmancmerge:prepmat'  ,  [ 'makmancmerge: input ' ,...
        'prep must be a file base name with no file type suffix' ]  )
      
    end % check prep
    
    % First input struct , if given , will be argument 2. The index before
    % that is 1.
    i = 1 ;
    
    % And build all necessary file names
    fdata = [  varargin{ 1 }  ,  '.mat'  ] ;
    fsort = [  varargin{ 1 }  ,  '.spksort.mat'  ] ;
    fout  = [  varargin{ 1 }  ,  '.manual.mat'  ] ;
    
  end % one or three inputs
  
  % Two or three inputs.
  if  any(  nargin  ==  [ 2 , 3 ]  )
    
    % There must be two structs , one after the other.
    if  ~ all ( cellfun(  @isstruct  ,  varargin( i + 1 : i + 2 )  ) )
    
      error (  'MAK:makmancmerge:ds'  ,  [ 'makmancmerge: input args ' ,...
          'd and s must be structs' ]  )
      
    end % check for structs
    
    % Extract structs
    d = varargin{ i + 1 } ;
    s = varargin{ i + 2 } ;
    
  end % two or three inputs
  
  % If there is only one input , prep , then the named data files must
  % exist
  if  nargin == 1
    
    % Look for files
    if  ~ all ( cellfun(  @( f ) exist( f , 'file' ) ,  ...
        { fdata , fsort }  ) )
    
      error (  'MAK:makmancmerge:nofiles'  ,  [ 'makmancmerge: input ' ,...
          'files %s and %s must exist' ]  ,  fdata  ,  fsort  )
      
    end % look for files
    
    % Load data
    load (  fdata  ,  'd'  )
    load (  fsort  ,  's'  )
    
  end % check input files
  
  % If there is no output requested then fout must be set to a valid file
  % name that is not in use
  if  nargout == 0
    
    % Do we have a file base name?
    if  isempty( fout )
      
      error (  'MAK:makmancmerge:noprep'  ,  [ 'makmancmerge: can''t' , ...
          ' write output file because input arg prep not given' ]  )
      
    % Does the output file exist yet?
    elseif  exist( fout , 'file' )
    
      error (  'MAK:makmancmerge:outfile'  ,  [ 'makmancmerge: output' ,...
          ' files already exist' ]  ,  fout  )
      
    end
    
  end % check output file
  
  % Redundant checking , if d or s is empty then there was a coding error
  if  isempty (  d  )  ||  isempty (  s  )
    
    error (  'MAK:makmancmerge:inputchk'  ,  [ 'makmancmerge: ' ,...
          'failed to get data structs d and s' ]  )
    
  end % check coding error
  
  
  %%% Prepare for manual merging %%%
  
  % Get master data structure
  m = getnewdata (  d  ,  s  ) ;
  
  % Make a merging tool
  h = makmergetool ;
  
  
  %%% Manual spike cluster merging %%%
  
  % Merge clusters on each electrode
  for  i = 1 : d.numtrodes
    
    % Electrode index
    eid = d.electrodes( i ) ;
    
    % Title header
    thead = sprintf (  'Electrode %d ( %d of %d )'  ,  ...
      eid  ,  i  ,  d.numtrodes  ) ;
    
    % Number of samples in a waveform
    ns = size (  s.wave_aligned{ i }  ,  1  ) ;
    
    % Get spike times for all spikes from this electrode , in PTB time
    t = cellfun (  @( s , t , e ) s + double( t( e  ==  eid ) )  ,...
      num2cell( d.start )  ,  d.spike.time  ,  d.spike.electrode  ,  ...
        'UniformOutput'  ,  false  ) ;
    
    % Collapse into a single vector
    t = [  t{ : }  ] ;
    
    % Reset man to empty string so that while loop performs at least one
    % iteration
    man = 'reset' ;
    
    % Reset loop runs so long as user hits a reset button
    while  ischar (  man  )
      
      % Null cutoff value
      jcut = -1 ;
    
      % Pass in data according to reset request
      switch  man
        
        % Pass in spike cluster assignments resulting from automated
        % merging
        case  'reset'
          
          Eargin = s.E{ i } ;
          spcarg = s.spikes_per_cluster{ i } ;
          cargin = s.clustind{ i } ;
            
        % Pass in initial spike cluster assignments prior to automated
        % merging
        case  'init'
          
          Eargin = s.E_init{ i } ;
          spcarg = s.spikes_per_cluster_init{ i } ;
          cargin = s.clustind_init{ i } ;
              
        % Redo automated merging up to a new connection strength cutoff
        % value
        case  'cutoff'
          
          % Get new cutoff value , recall that the previous iteration
          % returned man = 'cutoff' and E = value
          jcut = E ;
          
          % Redo merger with new cutoff value
          [ spcarg , cargin , Eargin ] = makcmerge (  ...
            s.spikes_per_cluster_init{ i }  ,  s.clustind_init{ i }  ,  ...
              s.E_init{ i }  ,  jcut  ) ;
      
      end % pass data
      
      % Perform manual mergers
      [ man , E , c , u ] = h.merge(  t  ,  s.wave_aligned{ i }  ,  ...
        s.pca{ i }  ,  Eargin  ,  spcarg  ,  cargin  ,  thead  ) ;
      
    end % reset loop
      
    % makmergetool closed
    if  all(  cellfun( @isempty , { man , E , c , u } )  )
      
      % Make sure that object is deleted
      delete (  h  )
      
      % Report and quit
      varargout{ 1 } = [] ;
      fprintf(  'MAK merge tool closed before all electrodes examined\n'  )
      return
      
    % New cutoff value used , known to be the case if jcut is not Null
    % value
    elseif  0  <=  jcut
      
      m.cutoff( i ) = jcut ;
      
    end % tool closed
    
    % List of mergers
    m.mergers{ i } = uint8 (  man  ) ;
    
    % Interface energy
    m.E{ i } = E ;
    
    % Number of clusters
    m.numclust( i ) = numel (  u  ) ;
    
    % Allocate cluster data structures
    cmap = zeros (  1  ,  size( E , 2 )  ,  'uint8'  ) ;
    wavg = zeros (  ns  ,  m.numclust( i )  ,  'single'  ) ;
    wvar = zeros (  ns  ,  m.numclust( i )  ,  'single'  ) ;
     rms = zeros (  1  ,  m.numclust( i )  ,  'single'  )  ;
   count = zeros (  m.numclust( i )  ,  d.numtrials  ,  'uint32'  ) ;
       K = false (  m.numclust( i )  ,  numel( c )  ) ;
    
    % Collect data for each cluster
    for  j = 1 : m.numclust( i )
      
      % Find spikes
      K( j , : ) = c  ==  u( j ) ;
      
      % Grab waveforms
      w = double (  s.wave_aligned{ i }( : , K( j , : ) )  ) ;
      
      % Compute average
      wavg( : , j ) = mean (  w  ,  2  ) ;
      
      % Compute variance
      wvar( : , j ) =  var (  w  ,  0  ,  2  ) ;
      
      % Reshape and ...
      w = reshape (  w  ,  numel( w )  ,  1  ) ;
      
      % Compute root-mean-squared
      rms( j ) = sqrt ( mean(  w .^ 2  ) ) ;
      
    end % clusters
    
    % Sort clusters based on rms ...
    [ rms , j ] = sort (  rms  ) ;
    
    % ... and apply sorting to all data
    u = u( j ) ;
    K = K( j , : ) ;
    wavg = wavg( : , j ) ;
    wvar = wvar( : , j ) ;
    
    % Build initial to final cluster index map
    cmap( u ) = 1 : m.numclust( i ) ;
    
    % Assign final cluster indices to spikes
    for  j = 1 : m.numclust( i )  ,  c( K(  j  ,  :  ) ) = j ;  end
    
    % Distribute spike cluster assignments to trials
    n = 0 ;
    
    for  j = 1 : d.numtrials
      
      % Find spikes on trial j
      k = n + 1 : n + d.spike.count( i , j ) ;
      n = n  +  d.spike.count( i , j ) ;
      
      % Find spikes from this electrode
      e = d.spike.electrode{ j }  ==  eid ;
      
      % Distribute spike cluster assignments
      m.clustind{ j }( e ) = c( k ) ;
      
      % Count spikes in each cluster
      count( : , j ) = arrayfun (  @( i ) sum( c( k )  ==  i )  ,  ...
        1 : m.numclust( i )  ) ;
      
    end % trials
    
    % Store data
    m.clustmap{ i } = cmap ;
    m.wave_avg{ i } = wavg ;
    m.wave_var{ i } = wvar ;
    m.rms{ i } = rms ;
    m.count{ i } = count ;
    
  end % electrodes
  
  
  %%% Final steps %%%
  
  % Close merge tool
  delete (  h  )
  
  % Save data
  if  ~ isempty (  fout  )
    
    % Report
    fprintf (  'Saving\n'  )
    
    % Write file
    save (  fout  ,  'm'  )
    
    % Remove write permissions , system command
    cstr = sprintf (  'chmod a-w %s'  ,  fout  ) ;
    
    % Execute write permission removal
    if  system (  cstr  ,  '-echo'  )
      
      error (  'MAK:makmancmerge:wperms'  ,  [ 'makmancmerge: ' , ...
        'failed to remove write permissions via system command:\n' , ...
          '  %s']  ,  cstr  )
      
    end % remove write permissions
    
  end % save data
  
  % Return data
  if  nargout  ,  varargout{ 1 } = m ;  end
  
  % Report
  fprintf (  'Done\n'  )
  
  
end % makmancmerge


%%% Subroutines %%%

% New data structure for etrode electrodes. Requires d and s input structs.
function  m = getnewdata (  d  ,  s  )
  

  % Number of electrodes
  nel = d.numtrodes ;
  
  
  % Ready list of connection strength cutoff values
  m.cutoff = s.J_cutoff ;
  
  % List of manual mergers per electrode
  m.mergers = cell (  nel  ,  1  ) ;
  
  % Final interface energy matrix following all mergers , initialised using
  % energy matrices that follow automated mergers
  m.E = s.E ;

  % The number of clusters that remain after all mergers
  m.numclust = zeros (  nel  ,  1  ,  'uint8'  ) ;
  
  % The mapping between the initial cluster indices from makprep and the
  % final cluster indices
  m.clustmap = cell (  nel  ,  1  ) ;
  
  % Average waveforms per cluster , per electrode
  m.wave_avg = cell (  nel  ,  1  ) ;
  
  % Waveform variance per cluster , per electrode
  m.wave_var = cell (  nel  ,  1  ) ;
  
  % Root-mean-squared per cluster , per electrode
  m.rms = cell (  nel  ,  1  ) ;
  
  % Spike count per cluster, trial, and electrode
  m.count = cell (  nel  ,  1  ) ;
  
  % Final cluster index of each spike on all electrodes per trial
  m.clustind = cellfun (  @( t ) zeros( size( t ) , 'uint8' )  ,  ...
    d.spike.time  ,  'UniformOutput'  ,  false  ) ;
  
  
end % getnewdata

