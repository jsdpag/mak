
function  [ R , TIEADJ ] = maktiedrank( X , varargin )
% 
% [ R , ... ] = maktiedrank( X , ... )
% 
% MET Analysis Kit. For some reason, Matlab's native tiedrank function only
% operates on vectors, not matrices. maktiedrank treats each column as X as
% a separate vector. So, tied ranks are found for each column of X,
% independently from all other columns. Accepts all the same inputs as
% tiedrank and returns all the same outputs. If X is a vector then
% maktiedrank acts identically to tiedrank.
% 
% Written by Jackson Smith - January 2020 - DPAG, University of Oxford
% 
  
  % X is a vector
  if  isvector( X )
    
    [ R , TIEADJ ] = tiedrank( X , varargin{ : } ) ;
  
    % Finished
    return
    
  end % X is a vector
  
  % Size of X in each dimension
  s = size( X ) ;
  
  % Number of columns
  n = prod( s( 2 : end ) ) ;
  
  % Allocate output, same type as X
  R = zeros( s , class( X ) ) ;
  
  % TIEADJ requested
  if  nargin > 1
    
    % Allocate TIEADJ
    TIEADJ = cell( 1 , n ) ;
    
    % Columns
    for  c = 1 : n
      [ R( : , c ) , TIEADJ{c} ] = tiedrank( X( : , c ) , varargin{ : } ) ;
    end
    
  % Just R is wanted
  else
    
    % Columns
    for  c = 1 : n
      R( : , c )  = tiedrank( X( : , c ) , varargin{ : } ) ;
    end
    
  end % cases
  
end % maktiedrank