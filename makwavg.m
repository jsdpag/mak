
function  a = makwavg (  x  ,  w  ,  varargin  )
% 
% a = makwavg (  x  ,  w  )
% a = makwavg (  x  ,  w  ,  dim  )
% a = makwavg (  __  ,  outtype  )
% a = makwavg (  __  ,  nanflag  )
% 
% MET Analysis Kit. Computes the average of data x weighted by w. x and w
% must be numeric arrays with the same size. The weighted average is
% computed as:
%   
%   a = [ w( 1 )x( 1 ) + w( 2 )x( 2 ) + w( 3 )x( 3 ) + ... ]  /  sum( w )
% 
% Optional input dim indicates along which dimension of x to compute the
% weighted average. outtype can be 'default', 'double', or 'native' while
% nanflag can be 'includenan' or 'omitnan'. See doc mean and doc sum for
% more information.
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% Input check %%%
  
  % Are x and w numeric?
  if  ~ isnumeric (  x  )  ||  ~ isnumeric (  w  )
    
    error (  'MAK:makwavg:xwnumeric'  ,  ...
      'makwavg: x and w must be numeric'  )
    
  % Are x and w the same size?
  elseif  numel (  x  )  ~=  numel (  w  )  ||  ...
          ndims (  x  )  ~=  ndims (  w  )  ||  ...
          any (  size(  x  )  ~=  size(  w  )  )
        
    error (  'MAK:makwavg:xwsize'  ,  ...
      'makwavg: x and w must have the same size'  )
    
  end % check input
  
  
  %%%  Weighted average %%%
  
  a = sum (  w  .*  x  ,  varargin{ : }  )  ./  ...
    sum (  w  ,  varargin{ : }  ) ;
  
  
end % makwavg

