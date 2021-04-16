
function  w = makpcawin (  n  ,  m  ,  p  )
% 
% w = makpcawin (  n  ,  m  ,  p  )
% 
% MET Analysis Kit, pre-processing. Returns a weighted window of n samples
% that is Gaussian-shaped, with a mean of m. The standard deviation of the
% Gaussian is calculated so that 100 * p percent of the area under the
% curve is reached by the nth sample. The window is normalised so as to sum
% to n. If p is zero then an n-element vector of ones is returned. w is a
% column vector. If n is zero then an empty matrix is returned.
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % p is zero , return ones instead. n is zero then return empty.
  if  ~ p
    w = ones ( n , 1 ) ;
    return
  elseif  ~ n
    w = [] ;
    return
  end
  
  % Find standard deviation for gaussian
  s = ( n  -  m )  /  norminv ( p ) ;

  % Get window's shape
  w = normpdf (  1 : n  ,  m  ,  s  ) ;

  % Normalise so that window sums to n
  w = w  /  sum ( w )  *  n ;
  
end % makpcawin

