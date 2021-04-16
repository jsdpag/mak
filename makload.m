
function  d = makload (  dname  ,  subid  ,  expid  ,  varargin  )
% 
% d = makload (  dname  ,  subid  ,  expid  )
% d = makload (  ...  ,  auxdir1  ,  auxdir2  ,  ...  auxdirN  )
% d = makload (  ...  ,  '-groupspikes'  )
% 
% MET Analysis Toolkit. Loads pre-processed data from directory dname for
% the subject with identifier subid, both of which are strings. Data is
% loaded from each experiment listed in the vector of identifiers expid.
% For subject subid and experiment id expid( i ), it is expected that
% directory dname will contain a pre-processed event data file from makprep
% containing variable 'd', and the results of manual spike sorting in a
% *.manual.mat file from makmancmerge containing variable 'm'.
% 
% A subset of data from m is added to d.cluster. Fields kept from m include
% .numclust, .wave_avg, .wave_std, .rms, and .count. These fields contain
% the number of clusters per electrode, the average cluster waveform, the
% waveform standard deviation, the waveform root-mean-squared, and the
% number of spikes for each trial per cluster. The field m.clustind is
% added to d.spike rather than d.cluster, because it contains spike cluster
% assignments per trial that is in register with other data in d.spike.
% 
% An optional list of auxiliary data directory names may follow the initial
% three input arguments. Each must be a string naming a directory that
% contains a file for each experiment that is named after the experiment's
% event data file in dname. The variables that these files contain are
% added to d.aux.(variable name). This is useful to add data that is highly
% task-specific.
% 
% The option flag '-groupspikes' may be given after all other input
% arguments. It instructs the function to group all spikes by spike
% cluster , using makspk. For the ith experiment , the default d( i ).spike
% data from makprep and makload is replaced with grouped data from makspk.
% Thus , d( i ).spike turns from a struct with fields: total , count ,
% time , electrode , channel , probe , unit , clustind into a struct with
% fields: time , elec , chan , probe , clus.
% 
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  
  %%% CONSTANTS %%%
  
  % Grouping option flag
  GRPFLG = '-groupspikes' ;
  
  % Fields to remove from manual spike sorting struct
  RMFNAM = { 'cutoff' , 'mergers' , 'E' , 'clustmap' , 'clustind' } ;
  
  
  %%% Input check %%%
  
  % Check that directory and subid arguments are all strings
  isstring (  dname  ,  'dname'  )
  isstring (  subid  ,  'subid'  )
  
  for  i = 4 : nargin  ,  isstring (  varargin{ i - 3 }  ,  i  )  ,  end
  
  % Experiment identifiers are all positive integers
  if  ~ isvector(  expid  )  ||  any(  expid < 1  |  ...
        ~ isreal( expid )  |  mod(  expid  ,  1  )  )
    
    error (  'MAK:makload:expid'  ,  ...
    'makload: input argument expid must be vector of positive integers'  )
    
  end % expid
  
  % Grouping option flag given
  if  3  <  nargin  &&  strcmp (  varargin{ end }  ,  GRPFLG  )
    
    % Raise flag
    GRPFLG = true ;
    
  % Flag not given
  else
    
    % Lower flag
    GRPFLG = false ;
    
  end % grouping flag
  
  
  %%% Search for file names %%%
  
  % Event data files
  dfile = cell ( size(  expid  ) ) ;
  
  % Manual spike sorting data files
  mfile = cell ( size(  expid  ) ) ;
  
  % Find pre-processed file names for each experiment
  for  i = 1 : numel (  expid  )
    
    % Search string for manual spike sorting file
    F = sprintf (  '%s.%d*.manual.mat'  ,  subid  ,  expid( i )  ) ;
    
    % Find experiment's manual spike sorting file
    F = dir ( fullfile(  dname  ,  F  ) ) ;
    
    if  isempty (  F  )
      
      ferr (  'manual' ,  'no manual spike sorting file' ,  subid ,  ...
        expid( i ) ,  dname  )
      
    elseif  1  <  numel (  F  )
      
      ferr (  'manual' ,  'too many manual spike sorting files' ,  ...
        subid ,  expid( i ) ,  dname  )
      
    end
    
    % Save manual spike sorting file name
    mfile{ i } = F.name ;
    
    % Get event file name
    dfile{ i } = strrep (  F.name  ,  '.manual'  ,  ''  ) ;
    
    if  ~ exist (  fullfile(  dname  ,  dfile{ i }  )  ,  'file'  )
      ferr (  'event' ,  'no event file' ,  subid ,  expid( i ) ,  dname  )
    end
    
    % Look for auxiliary files
    for  j = 1 : nargin - 3 - GRPFLG
      
      % File found , search for next
      if  exist (  fullfile( varargin{ j } , dfile{ i } )  ,  'file'  )
        continue
      end
      
      % File not found error
      ferr (  'aux' ,  'no auxiliary file' ,  subid ,  expid( i ) ,  ...
        varargin{ j }  )
      
    end % aux files
    
  end % file names
  
  
  %%% Load data %%%
  
  % Load each data struct separately then append
  d = cell ( size(  expid  ) ) ;
  
  % Experiments
  for  i = 1 : numel (  expid  )
    
    % Report
    fprintf (  'Loading  %s\n'  ,  strrep( dfile{ i } , '.mat' , '' )  )
    
    % Load event data
    F = fullfile (  dname  ,  dfile{ i }  ) ;
    d{ i } = load (  F  ,  'd'  ) ;
    d{ i } = d{ i }.d ;
    
    % Load manual spike sorting data
    F = fullfile (  dname  ,  mfile{ i }  ) ;
    m = load (  F  ,  'm'  ) ;
    m = m.m ;
    
    % Append manual spike sorting data to event data
    d{ i }.spike.clustind = m.clustind ;
    
    % Remove unwanted data from m
    m = rmfield (  m  ,  RMFNAM  ) ;
    
    % Store remaining cluster information
    d{ i }.cluster = m ;
    
    % There is auxiliary data , add field to event data struct
    if  3  <  nargin - GRPFLG  ,  d{ i }.aux = struct ;  end
    
    % Auxiliary data
    for  j = 1 : nargin - 3 - GRPFLG
      
      % Aux file name
      F = fullfile (  varargin{ j }  ,  dfile{ i }  ) ;
      
      % Load data
      a = load (  F  ) ;
      
      % Variable names
      avar = fieldnames (  a  ) ;
      
      if  isempty (  avar  )
        error (  'MAK:makload:auxdat'  ,  [ 'makload: auxiliary data ' ,...
          'file is empty %s' ]  ,  F  )
      end
      
      % Check whether variable names are in use
      k = isfield (  d{ i }.aux  ,  avar  ) ;
      
      if  any (  k  )
        error (  'MAK:makload:overloaded'  ,  [ 'makload: auxiliary ' ,...
          'variable names from %s are already in use: %s' ]  ,  F  ,  ...
            strjoin(  avar  ,  ' , '  )  )
      end
      
      % Add auxiliary data
      for  k = 1 : numel (  avar  )
        d{ i }.aux.( avar{ k } ) = a.( avar{ k } ) ;
      end
      
    end % aux datat
    
    % Do not group spikes by spike cluster , use default format. Continue
    % to next experiment.
    if  ~ GRPFLG  ,  continue  ,  end
    
    % Replace default .spike with grouped data
    d{ i }.spike = makspk (  d{ i }  ) ;
    
  end % experiments
  
  
  %%% Done %%%
  
  % Collapse into a struct vector
  d = [  d{ : }  ] ;
  
  
end % makload


%%% Subroutines %%%

% Check if argument is string , raise error if not
function  isstring (  s  ,  n  )
  
  % Argument is a string
  if  isrow(  s  )  &&  ischar(  s  )  ,  return  ,  end
  
  % If n is number then convert to string
  if  isnumeric (  n  )  ,  n = sprintf (  '%d'  ,  n  ) ;  end
  
  % Error , argument is not a string
  error (  'MAK:makload:notstring'  ,  ...
    'makload: input argument %s is not a string'  ,  n  )
  
end % isstring


% Raise an error indicating a missing data file
function  ferr (  eid  ,  s  ,  subid  ,  expid  ,  dname  )
  
  % Experiment id string
  eidstr = sprintf (  'MAK:makload:%s'  ,  eid  ) ;
  
  % Error message
  error (  eidstr  ,  [ 'makload: %s for subject %s, experiment %d, ' , ...
    'in directory %s' ]  ,  s  ,  subid  ,  expid  ,  dname  )
  
end % ferr

