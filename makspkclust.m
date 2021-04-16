
function  [ n , c , d0 ] = makspkclust (  par  ,  S  )
% 
% [ n , id , d0 ] = makspkclust (  par  ,  S  )
% 
% MAK Analysis Kit, pre-processing. Performs initial spike clustering on
% spike components in S. Returns vector n where n( i ) is the number of
% spikes in the ith cluster, and vector c where c( j ) is the cluster
% assignment of the jth spike. d0 is a scalar value ; it is the scaling
% term used to later compute interface energy.
% 
% Uses the initial clustering method of Fee et al. ( 1996 ; Hill et al.
% 2011 ). In brief, a large number of clusters are created relative to the
% number of neurones on the electrode. This over-clustering is useful
% because small clusters will track waveforms that slowly change shape, and
% yet they will be recombined according to powerful similarity measures.
% 
% Clusters are created following this outline:
% 
%   1) Sample 2xN cluster centres from existing N cluster centers by adding
%   two randomly sampled vectors to each centre
%   2) Assign each spike to the cluster with the nearest center.
%   3) Recompute each cluster's centre from its constituent spikes.
%   4) Repeat steps 2 & 3 several times so that cluster centres settle
%   down.
%   5) Repeat steps 1 to 4 until the desired number of bisections have
%   been performed.
% 
% It is advise to shuffle the random number generator's seed prior to
% calling this function for the first time.
% 
% 
% Input
% 
%   par - Struct parameter with fields:
%   
%     .bisecs - The number of cluster bisections to perform. Max of 7
%       produces up to 128 clusters.
%     .assign - The maximum number of spike cluster assignments to perform
%       before continuing to the next bisection.
%     .minspk - The minimum number of spikes required in a cluster. If a
%       cluster has too few spikes then those spikes are assigned to the
%       next-nearest cluster. This process repeats until all remaining
%       clusters have at least minspk spikes.
% 
%   S - N x M matrix of M spikes with N components each.
% 
% 
% Output
% 
% n - C element uint32 row vector with the number of spikes assigned to
%   each of C clusters.
% 
% c - M element uint8 row vector with the cluster assigment for each spike.
% 
% d0 - Scaling term for interface energy computation , single precision
%   floating point scalar value
% 
% 
% References:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% UltraMegaSort2000, https://neurophysics.ucsd.edu/software.php
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Input check %%%
  
  % Number of components per spike waveform and number of spikes
  [ ncmp , nspk ] = size (  S  ) ;
  
  % Maximum number of clusters
  cmax = 2 ^ par.bisecs ;
  
  % We cannot have fewer spikes than the minimum required by a cluster
  if  nspk  <  par.minspk
    
    error (  'MAK:makspkclust:S'  ,  ...
      'makspkclust: S must have at least %d columns'  ,  par.minspk  )
    
  % More clusters will be produced than can be indexed
  elseif  intmax ( 'uint8' )  <  cmax
    
    error (  'MAK:makspkclust:bisecs'  ,  [ 'makspkclust: up to %d ' , ...
      'clusters from %d bisections but max %d clusters supported' ]  ,  ...
      cmax  ,  par.bisecs  ,  intmax( 'uint8' )  )
    
  end % too few spikes
  
  
  %%% Preparation %%%
  
  % In order to find the scale of randomness to add to new cluster centres,
  % estimate the average distance between spikes by sampling 5000 pairs
  Crnd = estdist (  5e3  ,  S  ) ;
  
  % Heuristic randomness scaling factor used by UltraMegaSort2000
  Crnd = Crnd  /  100  /  ncmp ;
  
  % Allocate output vectors
  n = zeros (  1  ,  cmax  ,  'uint32'  ) ;
  c = zeros (  1  ,  nspk  ,  'uint8'   ) ;
  
  % Store old cluster assignment to see if these have changed
  oldc = zeros (  1  ,  nspk  ,  'uint8'   ) ;
  
  % Cluster centres
  cen = zeros (  ncmp  ,  cmax  ,  'single'  ) ;
  
  % Initialise first cluster
  cnum = 1 ;
  n( 1 ) = nspk ;
  cen( : , 1 ) = mean (  S  ,  2  ) ;
  ci = 1 : cnum ;
  
  
  %%% Clustering %%%
  
  % Bisections
  for  bisecs = 1 : par.bisecs
    
    
    %- New clusters -%
    
    % There is no old cluster assignment , guarantee that assignment loop
    % won't break by accidental equality of new c values and oldc from the
    % previous bisection iteration
    oldc( : ) = 0 ;
    
    % Update number of cluster centres
    cnum = 2  *  numel (  ci  ) ;
    
    % Duplicate existing cluster centres
    cen( : , 1 : cnum ) = repmat (  cen( : , ci )  ,  1  ,  2  ) ;
    
    % New cluster index
    ci = 1 : cnum ;
    
    % Add random noise to each cluster centre
    cen( : , ci ) = Crnd * rand (  ncmp  ,  cnum  )  +  cen( : , ci ) ;
    
    
    %- Assign spikes to clusters -%
    
    % Do this a maximum number of times
    for  assign = 1 : par.assign
      
      % Compute distances of each spike to each cluster centre.
      % cen( : , ci ) can return fewer than cnum centres if clusters were
      % dropped in a previous assignment iteration due to lack of spikes.
      d = pdist2 (  S'  ,  cen( : , ci )'  ) ;
      
      % Assign spikes to nearest clusters ...
      [ ~ , c( : ) ] = min (  d  ,  []  ,  2  ) ;
      
      % ... and readjust from d column index to cluster index
      c( : ) = ci (  c  ) ;
      
      % The number of spikes per cluster
      n( ci ) = arrayfun (  @( i ) sum( c  ==  i )  ,  ci  ) ;
      
      % Look for clusters with too few spikes
      for  i = ci(  n( ci )  <  par.minspk  )
        
        % It is possible that this cluster now has enough spikes after
        % re-assigning those from another cluster that was too small. Go to
        % the next cluster in this case.
        if  par.minspk  <=  n( i )  ,  continue  ,  end
        
        % Find all spikes belonging to this cluster
        s = find (  c  ==  i  ) ;
        
        % We will need to search for the next closest cluster. To guarantee
        % that no spike is re-assigned to cluster i, we make all the
        % spike distances to this cluster infinite.
        d( : , ci == i ) = Inf ;
        
        % Now we can look for the next closest clusters
        [ ~ , c( s ) ] = min (  d( s , : )  ,  []  ,  2  ) ;
        
        % Remember to readjust column to cluster indices
        c( s ) = ci( c( s ) ) ;
        
        % No spikes remain in the current cluster
        n( i ) = 0 ;
        
        % Count the number of spikes assigned to each new cluster
        for  j = unique( c( s ) )
          n( j ) = sum( c( s )  ==  j )  +  n( j ) ;
        end
        
      end % too few spikes
      
      % Find empty clusters and remove them from index vector
      ci(  n( ci )  ==  0  ) = [] ;
      
      % No change in cluster assignment , stop assigning
      if  all (  oldc  ==  c  )  ,  break  ,  end
      
      % Remember new cluster assignments to compare against future
      % assignments
      oldc( : ) = c ;
      
      % Re-calculate cluster centres
      for  i = ci
        
        % Find spikes in this cluster
        j = c  ==  i ;
        
        % New cluster centre
        cen( : , i ) = mean (  S( : , j )  ,  2  ) ;
        
      end % new centres
      
    end % assign spks
    
  end % bisections
  
  
  %%% Finalise output %%%
  
  % Discard empty clusters
  n = n( ci ) ;
  cen = cen( : , ci ) ;
  
  % If there were empty clusters ...
  if  numel (  ci  )  <  cnum
  
    % ... then re-label clusters and spike assignments
    for  i = 1 : numel ( ci )

      % Find all spikes in cluster
      j = c  ==  ci( i ) ;

      % Re-label
      c( j ) = i ;

    end % re-label
    
  end % zero-spike clusters
  
  % Scaling term used by UltraMegaSort2000 , line one is W = T - B , line
  % two acts on W
  d0 = cov( S' ) - cov( cen( : , c )' ) ;
  d0 = sqrt ( sum( diag(  d0  ) ) )  /  10 ;
  
  
end % makspkclust


%%% Sub-routines %%%

% Estimate the mean distance between pairs of spikes by sampling n pairs
% of different spikes from S
function  mdist = estdist (  n  ,  S  )
  
  % Build up a distribution of pairwise distances
  d = zeros (  n  ,  1  ) ;
  
  % Count the number of distances found between different spikes
  ndist = 0 ;
  
  % Number of spikes
  nspk = size (  S  ,  2  ) ;
  
  % Sample distances until the desired number are found
  while  ndist  <  n
    
    % Randomly sample column indices for S
    i = ceil (  nspk  *  rand (  2  ,  n - ndist  )  ) ;
    
    % Find pairs that compare a spike with itself
    j = i( 1 , : )  ==  i( 2 , : ) ;
    
    % Remove them
    i(  :  ,  j  ) = [] ;
    
    % Number of new distances
    nnew = size (  i  ,  2  ) ;
    
    % All pairs removed , try again
    if  ~ nnew  ,  continue  ,  end
    
    % Index vector for distance buffer where new distances will go
    k = ndist + 1  :  ndist + nnew ;
    
    % Calculate distances
    d( k ) = euclid (   S( : , i( 1 , : ) ) ,   S( : , i( 2 , : ) )   ) ;
    
    % Update number of distances found
    ndist = ndist  +  nnew ;
    
  end % sample distances
  
  % Mean distance between spikes
  mdist = mean (  d  ) ;
  
end % estdist


% Euclidean distance
function  d = euclid (  s  ,  t  )
  
  d = sqrt ( sum(  (  s  -  t  ) .^ 2  ,  1  ) ) ;
  
end % euclid

