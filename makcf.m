
function  varargout = makcf( varargin )
% 
% [ ... ] = makcf( ... )
% 
% MET Analysis Kit. Functions the same way as cellfun, but 'UniformOutput'
% is always false. Thus, the output is always a cell array of objects. In
% other words, this is the same as running:
% 
%   [ ... ] = cellfun( ... , 'UniformOutput' , false )
% 
% Written by Jackson Smith - January 2020 - DPAG, University of Oxford
% 
  
  % Allocate output
  varargout = cell( 1 , nargout ) ;
  
  % Evaluate
  [ varargout{ : } ] = cellfun( varargin{ : } , 'UniformOutput' , false ) ;
  
end % makcf
