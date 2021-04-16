
function  r = makgabor (  c  ,  x  )
% 
% r = makgabor (  c  ,  x  )
% c = makgabor
% 
% MET Analysis Kit. Evaluates a 1-dimensional Gabor function with
% coefficients at values in x. Negative values are half-wave rectified to
% zero. c can be either a six-element double vector of real values or a
% struct in which all fields contain a scalar real double ; the struct must
% have fields:
%   
%   y0 - Baseline value
%   A  - Gaussian envelope amplitude
%   x0 - Gaussian envelope horizontal offset
%   s  - Width of the Gaussian envelope (s for sigma)
%   f  - Frequency of cosine
%   p  - Phase of cosine , in radians
% 
% If c is a vector then it must contain the same coefficients in the
% order: [ y0 , A , x0 , s , f , p ]. The vector form makes makgabor
% suitable for use with fitting functions such as lsqcurvefit.
% 
% If no inputs are provided then a struct version of c is returned with
% fields initialised to zero.
% 
% For an example equation and procedure for fitting disparity tuning curves
% see:
% 
%   Tanabe S, Umeda K, Fujita I. Rejection of false matches for binocular
%     correspondence in macaque visual cortical area V4. J Neurosci. 2004
%       Sep 15;24(37):8170-80.
% 
% Written by Jackson Smith - April 2018 - DPAG , University of Oxford
% 
  
  
  %%% No input %%%
  
  % Make c a vector of zeros that will be used below
  if  ~ nargin  ,  c = zeros (  1  ,  6  ) ;  end
  
  
  %%% Coefficients %%%
  
  % Coefficient vector given , turn this into a struct
  if  ~ isstruct (  c  )  &&  isvector (  c  )
    
    c = struct (  'y0' ,  c( 1 ) ,  'A' ,  c( 2 ) ,  'x0' ,  c( 3 ) ,  ...
      's' ,  c( 4 ) ,  'f' ,  c( 5 ) ,  'p' ,  c( 6 )  ) ;
    
    % No input , return c and quit
    if  ~ nargin
      r = c ;
      return
    end
    
  end % vect to struct
  
  
  %%% Evaluate gabor %%%
  
  % Data minus Gaussian offset
  dx = x  -  c.x0 ;
  
  % Sigma from standard deviation to twice the variance
  c.s = 2  *  c.s ^ 2 ;
  
  % Get the Gaussian enveloped cosine
  r = exp(  - dx .^ 2  /  c.s  )  .*  cos(  2 * pi * c.f * dx  +  c.p  ) ;
  
  % Scale and add baseline
  r = c.A * r  +  c.y0 ;
  
  % And the good lord said: Rectify, my Brothers and Sisters. Rectify!
  r = max (  r  ,  0  ) ;
  
  
end % makgabor

