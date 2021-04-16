
function  [ s , e , d ] = makskiptime (  t  ,  m  )
% 
% [ s , e , d ] = makskiptime (  t  ,  m  )
% 
% MET Analysis Kit. Computes the start ( s ) and end time ( e ) of each
% skipped i.e. missed frame, and the estimated duration ( d ). s and e are
% in seconds from the start of the trial, and d is in seconds. Frame onset
% times are given in vector t, in seconds from the start of the trial. m is
% a vector of equal size to t. m( i ) is non-zero if t( i ) is the onset
% time of a skipped frame , meaning that the frame buffer flip occurred
% after a given deadline. If the ith frame is the jth skipped frame, then
% s( j ), e( j ), and d( j ) are computed from frame times t( i - 1 ) and
% t( i ). Returns empties if there were no skipped frames.
% 
% Written by Jackson Smith - June 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % t and m are vectors of equal size
  if  ~ isvector (  t  )  ||  ~ isvector (  m  )  ||  ...
      numel (  t  )  ~=  numel (  m  )
    
    error (  'MAK:makskiptime:tmvector'  ,  ...
      'makskiptime: t and m must be vectors of equal size'  )
    
  % t must be numeric, real-valued, and finite
  elseif  ~ isnumeric (  t  )  ||  ~ isreal (  t  )  ||  ...
      ~ all ( isfinite(  t  ) )
    
    error (  'MAK:makskiptime:tnum'  ,  ...
      'makskiptime: t must be a real-valued, finite, numeric vector'  )
    
  % m must be numeric, real-valued, and finite or logical
  elseif  ~ ( isnumeric (  m  )  ||  islogical (  m  ) )  ||  ...
      ~ isreal (  m  )  ||  ~ all ( isfinite(  m  ) )
    
    error (  'MAK:makskiptime:mnum'  ,  [ 'makskiptime: m must be a ' , ...
      'real-valued, finite, numeric or logical vector' ]  )
    
  end % check input
  
  % Guarantee that m is logical
  if  ~ islogical (  m  )  ,  m = m  ~=  0 ;  end
  
  % Guarantee that t and m are column vectors
  if  ~ iscolumn (  t  )  ,  t = t( : ) ;  end
  if  ~ iscolumn (  m  )  ,  m = m( : ) ;  end
  
  
  %%% Compute skip onset and end times %%%
  
  % Default return values
  s = [] ;  e = [] ;  d = [] ;
  
  % Find skipped frames
  i = find (  m  ) ;
  
  % No skipped frames
  if  isempty (  i  )  ,  return  ,  end
  
  % Look for special case when first frame was skipped
  if  i( 1 )  ==  1
    
    % Start time is zero and end time is estimated onset time of that frame
    s = 0 ;  e = t( 1 ) ;
    
    % This frame is accounted for , drop it from the list
    i( 1 ) = [] ;
    
  end % special case
  
  % Start times
  s = [  s  ;  t( i - 1 )  ] ;
  
  % End times
  e = [  e  ;  t( i )  ] ;
  
  % Durations
  d = e  -  s ;
  
  
end % makskiptime

