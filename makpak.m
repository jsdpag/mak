
function  y = makpak (  r  ,  h  ,  varargin  )
%
% y = makpak (  r  ,  h  ,  x1  ,  x2  ,  ...  xN  )
% 
% MET Analysis Kit. Returns M output arguments from function h in y , an M
% by 1 element cell array. Input arguments x1 to xN are passed directly to
% function h. r is a vector of indices naming which output arguments from h
% to return , in the order specified. If r is a scalar then the output from
% h is returned directly as y , without wrapping in a cell array.
% 
% makpak is indended for use with makfun when outputs other than return
% argument 1 are required. The reason that makfun cannot return multiple
% outputs is that it can operate on multiple functions that may have
% different limits on their maximum number of output arguments.
% 
% Example: Return the rccg and nscx output arguments from makrccg , for the
%   first 0.5 second of responses of each trial with trials grouped by
%   baseline disparity , for each spike cluster. And return nscx first.
% 
%   Let R be the T by C cell array of spike times over T trials and C spike
%   clusters. Let D be a T-element vector of baseline disparity values on
%   each trial.
%   
%   % Make function handle that packs makrccg output arguments into a 
%   % single cell array
%   h = @( r ) makpak( [ 3 , 1 ] , makrccg , [ 0.044 , 0.544 ] , r ) ;
%   
%   % Return multiple makrccg outputs for each trial grouping
%   [ y , n , u ] = makfun (  h  ,  R  ,  D  ) ;
%   
%   % Obtain each input argument over all trial groups
%    ccg = cell2mat (  y( 1 , : )  ) ;
%   rccg = cell2mat (  y( 2 , : )  ) ;
%
%   If M = numel (  n  ) then y is a 2 by M cell array. Row 1 of y contains
%   nscx for each trial grouping and row 2 contains rccg.
%   
% Written by Jackson Smith - May 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Number of input and output arguments
   narginchk (  2  ,  Inf  )
  nargoutchk (  0  ,    1  )
  
  % r needs to be a numeric vector with real integer values
  if  ~ isvector (  r  )  ||  ~ isnumeric (  r  )  ||  ...
      ~ all( isfinite(  r  ) )  ||  ~ isreal (  r  )  ||  ...
      any (  r  <  1  |  mod(  r  ,  1  )  )
    
    error (  'MAK:makpak:r'  ,  [ 'makpak: r must be a finite ' , ...
      'numeric vector of real integer values greater than zero' ]  )
    
  % h must be a function handle
  elseif  ~ isa (  h  ,  'function_handle'  )
    
    error (  'MAK:makpak:h'  ,  'makpak: h must be a function handle'  )
    
  end % input checks
  
  
  %%% Return output arguments %%%
  
  % Maximum output argument index
  M = max (  r  ) ;
  
  % Allocate return cell array
  y = cell (  M  ,  1  ) ;
  
  % Execute function and get return values
  [  y{ : }  ] = h (  varargin{ : }  ) ;
  
  % Get named output arguments
  y = y( r ) ;
  
  % But there was only one , so return that directly
  if  isscalar (  y  )  ,  y = y{ 1 } ;  end
  
  
end % makpak

