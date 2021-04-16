
function  spk = makspk (  d  )
% 
% spk = makspk (  d  )
% 
% MET Analysis Kit. Groups all spike times by trial and by spike cluter ,
% for each experiment in d. d must be a struct array of the kind returned
% by makload. Returns spk , a struct with the same size as d. It has
% fields:
%   
%   .time - An N x C cell array. Rows are indexed over N trials and columns
%     are indexed over C spike clusters. Each element contains a single
%     floating point row vector of spike times for a specific cluster on a
%     specific trial. Times are in seconds from the start of the analysis
%     epoch. Column indexing progresses from the lowest to highest cluster
%     index on each electrode, and second from the lowest to highest
%     electrode ID. Thus, spike clusters from an electrode are ordered
%     consecutively, and sets of clusters are ordered by electrode. If an
%     electrode has no spike cluster then no column refers to it.
%   
%   .elec - uint8 row vector gives the electrode ID of each column of
%     .time.
%   
%   .chan - uint16 row vector gives the channel ID of each column of .time.
%   
%   .probe - uint8 row vector gives the probe ID of each column of .time.
%   
%   .clus - uint8 row vector gives the cluster ID of each column of .time.
%   
% Hence , to examine spike times from cluster i on electrode j , we get the
% column index of .time as k = spk.clus == i  &  spk.elec == j. Spikes for
% all trials from that specific cluster are retrieved by spk.time( : , k ),
% a cell array column vector.
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford    
% 
  
  
  %%% Constants %%%
  
  % Return cell arrays from *fun ( ) function calls. Pass in as comma-
  % separated list.
  UF = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Check input %%%
  
   narginchk (  1  ,  1  )
  nargoutchk (  0  ,  1  )
  
  if  ~ isstruct (  d  )
    error (  'MAK:makspk:inputd'  ,  'makspk: d must be a struct array'  )
  end
  
  
  %%% Allocate output struct %%%
  
  % We get a cell array of empties that is the same size as d
  spk = cell ( size(  d  ) ) ;
  
  % And use this as the value for each field to get a struct array the size
  % of d. Each field contains empty.
  spk = struct (  'time' ,  spk ,  'elec' ,  spk ,  'chan' ,  spk ,  ...
    'probe' ,  spk  ,  'clus' ,  spk  ) ;
  
  
  %%% Group spike times %%%
  
  % Experiments
  for  i = 1 : numel (  d  )
    
    % Electrode ID set
    eid = d( i ).electrodes ;
    
    % Number of clusters per electrode , as row vector
    nclust = d( i ).cluster.numclust' ;
    
    % Spike times per trial
    T = d( i ).spike.time ;
    
    % Electrode IDs per trial
    E = d( i ).spike.electrode ;
    
    % Cluster ID per trial
    C = d( i ).spike.clustind ;
    
    % Find electrodes with clusters
    j = 0  <  nclust ;
    
    % Get electrode ID per column
    spk( i ).elec = arrayfun (  @( e , n ) repmat( e , 1 , n )  ,  ...
      eid( j )  ,  nclust( j )  ,  UF{ : }  ) ;
    spk( i ).elec = cell2mat (  spk( i ).elec  ) ;
    
    % Map to channel and probe IDs
    spk( i ).chan  = d( i ).elec2chan (  spk( i ).elec  ) ;
    spk( i ).probe = d( i ).elec2probe(  spk( i ).elec  ) ;
    
    % Get cluster index of each column
    spk( i ).clus = arrayfun ( @( n ) 1 : n  ,  nclust  ,  UF{ : }  ) ;
    spk( i ).clus = cell2mat (  spk( i ).clus  ) ;
    
    % Allocate cell array
    spk( i ).time = cell (  numel( T )  ,  sum( double( nclust ) )  ) ;
    
    % Column index for .time
    ci = 0 ;
    
    % Electrodes
    for  j = 1 : d( i ).numtrodes
      
      % Spike clusters
      for  k = 1 : nclust( j )
        
        % Increment column index
        ci = ci  +  1 ;
        
        % Group spike times on all trials from this cluster
        spk( i ).time( : , ci ) = cellfun (  ...
          @( t , e , c ) t(  e == eid( j )  &  c == k  ) ,...
            T ,  E ,  C ,  UF{ : }  ) ;
        
      end % clusters
    end % electrodes
  end % experiments
  
  
end % makspk

