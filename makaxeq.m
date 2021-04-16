
function  lim = makaxeq (  varargin  )
% 
% lim = makaxeq
% lim = makaxeq (  ax  )
% 
% MET Analysis Kit. With no input args, takes the current axes and sets it
% to be square and tight around the data. The limits of the x- and y-axis
% are then matched to the minimum and maximum value observed along either
% one. This makes each axis directly comparable to the other in the plot,
% and the 45 degree line from the lower-left to upper-right corners will be
% the equality line.
% 
% ax can be a single handle to a specific axes, or an array of handles to
% different axes. The axis equalisation is applied to each given axes such
% that all of them use the same axis limits, for direct comparisons of data
% between different axes.
% 
% The axis limits applied in all cases is returned in lim.
% 
% Written by Jackson Smith - Sept 2019 - University of Oxford, DPAG
% 
  
  % Input/output check
  narginchk  ( 0 , 1 )
  nargoutchk ( 0 , 1 )
  
  % Input given?
  if  nargin
    
    % Fetch input
    ax = varargin{ 1 } ;
    
    % Check that we have graphics objects
    if  any ( ~ isgraphics( ax ) )
      
      error (  'Non-graphics object found in ax'  )
      
    end % graphics objects
    
    % Get axes
    ax = findobj (  ax  ,  'Type'  ,  'axes'  ) ;
    
    % Nothing left
    if  isempty (  ax  )
      
      error (  'No axes found in ax'  )
      
    end % no input
    
  else
    
    % No input, look for currect axes
    ax = gca ;
    
  end % input check
  
  % Set all axes to be square in shape and to match limits to the data
  axis (  ax  ,  'square'  ,  'tight'  )
  
  % x- and y-axis limits
  xlim = get (  ax  ,  'XLim'  ) ;
  ylim = get (  ax  ,  'YLim'  ) ;
  
  % Global limits
  if  isscalar (  ax  )
    lim = [ xlim , ylim ] ;
  else
    lim = [ xlim{ : } , ylim{ : } ] ;
  end
  
  lim = [ min( lim ) , max( lim ) ] ;
  
  % Set limits
  axis (  ax  ,  [ lim , lim ]  )
  
end % makaxeq

