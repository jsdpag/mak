
function  fig = makfig( varargin )
% 
% fig = makfig( mul )
% fig = makfig( fig )
% fig = makfig( fig , mul )
%
% MIND Analysis Kit. Returns a figure handle after changing the size of the
% figure to be mul * [ 21.0 , 29.7 ] centimetres. mul is optional, with
% a default value of 1, producing a figure the size of a piece of A4 paper.
% If mul < 0 then a landscape orientation is used instead of portrait.
% Optionally, fig can be an existing, valid figure handle. A new figure is
% created by default, and its handle is returned. If a handle is provided
% then the same one is returned, after changing the figure's size.
%
% It's worth pointing out that at best this will be an approximation.
% 
% Created by Jackson Smith - April 2021 - ESI (Fries Lab)
% 


%%% Constants %%%

% Size of A4 in centimetres, portrait orientation [ width , height ]
A4 = [ 21.0 , 29.7 ] ;

% Default multiplier
DEFMUL = 1 ;

% Point to Root graphics object
GRT = groot ;


%%% Check input %%%

% Number of inputs/outputs
 narginchk( 0 , 2 )
nargoutchk( 0 , 1 )

% No input
if  nargin == 0
  
  % Defaults
  fig = figure ;
  mul = DEFMUL ;
  
% Just one input
elseif  nargin == 1
  
  % This is a figure
  if  isgraphics( varargin{ 1 } , 'figure' )
    
    % Assign figure
    fig = varargin{ 1 } ;
    
    % Default multiplier
    mul = DEFMUL ;
    
  % Otherwise it must be the multiplier
  else
    
    % Assign multiplier
    mul = varargin{ 1 } ;
    
    % Create figure
    fig = figure ;
    
  end % assign input
  
% Two inputs
else
  
  % Assign names to input args
  [ fig , mul ] = varargin{ : } ;
  
end % check/get input

% Check figure handle and multiplier for correctness
if  ~ isscalar( fig ) || isnumeric( fig ) || ~isgraphics( fig , 'figure' )
  
   error( 'MAK:makfig:fig' , ...
     'makfig: fig must be scalar, valid figure handle' )
  
elseif  ~isscalar( mul ) || ~isnumeric( mul ) || ~isfinite( mul ) || ~mul
  
  error( 'MAK:makfig:fig' , ...
     'makfig: mul must be scalar, finite, non-zero number' )
  
end % check input is correct

% Guarantee mul is double
if  ~ isa( mul , 'double' ) , mul = double( mul ) ; end

% mul is less than zero
if  mul < 0
  
  % Reverse orientation of figure from portrait to landscape
  A4 = flip( A4 ) ;
  
  % Take absolute value of multiplier
  mul = - mul ;
  
end % portrait -> landscape


%%% Scale figure %%%

% Current unit type of figure
unit = fig.Units ;

% Switch to centimetres
fig.Units = 'centimeters' ;

% Change size of drawable area
fig.Position( 3 : 4 ) = mul .* A4 ;

% Now make units equal to root graphics object
fig.Units = GRT.Units ;

% Is the top of the figure visible?
if  sum( fig.OuterPosition( [ 2 , 4 ] ) )  >  GRT.ScreenSize( 4 )
  
  % No, move figure down until the top is visible
  fig.OuterPosition( 2 ) = GRT.ScreenSize( 4 ) - fig.OuterPosition( 4 ) ;
  
end % guarantee visibility of top

% Restore original unit
fig.Units = unit ;


%%% DONE - makfig %%%

