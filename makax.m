
function  ax = makax( varargin )
% 
% ax = makax( ... )
% ax = makax( ax , ... )
% 
% MIND Analysis Kit. Returns a new axes object handle in ax. This is
% automatically formatted with options that may be useful for generating
% publication quality plots.
% 
% Accepts any set of input arguments that is also accepted by the axes( )
% function (see doc axes for details). When input is provided then axes is
% called as follows:
% 
%   ax = axes( < USER INPUT > , < DEFAULT PARAMETERS > )
% 
% User input name/value pairs will overide default settings.
% 
% If first argument ax is an existing axes object then user-defined and
% default parameters are applied to it. Returns same axes handle.
% 
% Default axes parameters are:
% 
%      parent - Current figure
%     TickDir - 'out'
%     TickLen - [ 0.025 , 0.025 ]
%   LineWidth - 1.0
%      XColor - 'k'      Same for Y & Z
%    XLimSpec - 'tight'  Same for Y & Z
%    XLimMode - 'auto'   Same for Y & Z
%    NextPlot - 'add'
%         Box - 'off'
%    FontSize - 10
%   FontWeight - 'normal'
%     FontName - 'Helvetica'
% 
% Default axes label parameters (see ax.XLabel, etc.)
% 
%     FontSize - 11
%   FontWeight - 'normal'
%     FontName - 'Helvetica'
%        Color - 'k'
% 
% Created by Jackson Smith - April 2021 - ESI (Fries Lab)
% 


%%% CONSTANTS %%%

% Default axes parameters. Each column is a name/value pair, name in row 1
% and parameter value in row 2.
DEFPAR = {  'TickDir' , 'out' ;
            'TickLen' , [ 0.025 , 0.025 ] ;
          'LineWidth' , 1.0 ;
             'XColor' , 'k' ;
             'YColor' , 'k' ;
             'ZColor' , 'k' ;
           'XLimSpec' , 'tight' ;
           'XLimMode' , 'auto' ;
           'YLimSpec' , 'tight' ;
           'YLimMode' , 'auto' ;
           'ZLimSpec' , 'tight' ;
           'ZLimMode' , 'auto' ;
           'NextPlot' , 'add' ;
                'Box' , 'off' ;
           'FontSize' , 10 ;
         'FontWeight' , 'normal' ;
           'FontName' , 'Helvetica' }' ;

% Axis label parameters
AXIPAR = { 'FontSize' , 11 , 'FontWeight' , 'normal' , ...
  'FontName' , 'Helvetica' , 'Color' , 'k' } ;


%%% Check number of outputs %%%

nargoutchk( 0 , 1 )


%%% Get axes handle %%%

% User has provided at least one input. First item in the list is a valid
% axes object handle.
if  nargin  &&  isgraphics( varargin{ 1 } , 'axes' )
  
  % Point to axes
  ax = varargin{ 1 } ;
  
  % Apply user input arguments
  if  nargin > 1 , set( ax , varargin{ 2 : end } ) ; end
  
% Axes not provided
else
  
  % Create axes object, and apply user input args
  ax = axes( varargin{ : } ) ;
  
end % get axes handle


%%% Apply user input first %%%

% User input arguments containing character strings
for  I = varargin , i = I{ 1 } ;
  
  % Input is not a string, skip to next
  if  ~ ischar( i )  ||  ~ isrow( i ) , continue , end
  
  % See if this names a default axes parameter
  j = strcmp( DEFPAR( 1 , : ) , i ) ;
  
  % This is not a default axes parameter, continue
  if  ~ any( j ) , continue , end
  
  % User input overides default value
  DEFPAR( : , j ) = [ ] ;
  
end % input args


%%% Apply default parameters %%%

% Apply default axes parameters
set( ax , DEFPAR{ : } ) ;

% Axis labels
for  I = { 'Title' , 'XLabel' , 'YLabel' , 'ZLabel' } , i = I{ 1 } ;
  
  % Apply default label parameters
  set( ax.( i ) , AXIPAR{ : } ) ;
  
end % labels


%%% DONE - makax %%%

