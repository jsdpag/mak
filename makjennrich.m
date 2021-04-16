
function  [ pval , chi2 , df ] = makjennrich( R1 , n1 , R2 , n2 )
% 
% [ pval , chi2 , df ] = makjennrich( R1 , n1 , R2 , n2 )
% 
% MET Analysis Kit. An implementation of the Jennrich Test for the equality
% of two correlation matrices. The test is applied to correlation matrices
% R1 and R2, each computed from n1 and n2 samples, respectively. Samples
% for both matrices are assumed to be independent and taken from p-variate
% normal distributions. Use another test e.g. Lilliefors test (lillietest)
% to check if the data used to compute R1 and R2 have normal distributions.
% 
% Returns the p-value of the test (pval), the chi-squared value (chi2), and
% the degrees of freedom (df).
% 
% Variable names conform, within the limits of the ASCII character set, to
% conventions used by Jennrich (1970). Input variables must all be of a
% numeric type.
% 
% Reference
% 
%   Jennrich, R. I. (1970). "An Asymptotic chi-square Test for the Equality
%     of Two Correlation Matrices." Journal of the American Statistical
%     Association 65(330): 904-912.
% 
% See also
% 
%    Maximilian A.M. VERMORKEN (2020). Jennrich Test (https://
%      www.mathworks.com/matlabcentral/fileexchange/15171-jennrich-test),
%      MATLAB Central File Exchange. Retrieved December 18, 2020.
% 
% Written by Jackson Smith - Décembre 2020 - DPAG , University of Oxford
% 
  
  %%% Check Input %%%
  
  % Collect number of dimensions (number of cols/rows in square matrix)
  p = [ 0 , 0 ] ;
  
  % Correlation matrices
  for  R = { { R1 , n1 , 1 } , { R2 , n2 , 2 } }
    
    % Give meaningful names 'm'atrix, 'n'umber of samples, 'i'ndex number
    [ m , n , i ] = R{ 1 }{ : } ;
    
    % No empty data allowed
    if  isempty( m )
      
      error( 'MAK:makjennrich:empty_matrix' , ...
        'makjennrich: R%d is empty' , i )
      
    elseif  isempty( n )
      
      error( 'MAK:makjennrich:empty_samples' , ...
        'makjennrich: n%d is empty' , i )
      
    end
    
    % Is this a numeric array?
    if  ~ isnumeric( m )
      
      error( 'MAK:makjennrich:numeric_matrix' , ...
        'makjennrich: R%d is not a numeric type' , i )
      
    % Is this a matrix?
    elseif  ~ ismatrix( m )
      
      error( 'MAK:makjennrich:matrix' , ...
        'makjennrich: R%d is not a 2D matrix' , i )
      
    % Are all elements a valid value?
    elseif  ~ all( isfinite( m( : ) ) )
      
      error( 'MAK:makjennrich:numeric' , ...
        'makjennrich: R%d contains invalid values i.e. not finite' , i )
      
    end % check matrix
    
    % Get number of rows and columns in the matrix
    [ rows , cols ] = size( m ) ;
    
    % Not a square matrix
    if  rows ~= cols
      
      error( 'MAK:makjennrich:square' , ...
        'makjennrich: R%d is not a square matrix' , i )
      
    end % check matrix
    
    % Store dimensionality of data used to compute Ri
    p( i ) = rows ;
    
    % Is number of samples scalar
    if  ~ isscalar( n )
      
      error( 'MAK:makjennrich:scalar' , ...
        'makjennrich: n%d is not scalar' , i )
      
    % Is number of samples a numeric type?
    elseif  ~ isnumeric( n )
      
      error( 'MAK:makjennrich:numeric_samples' , ...
        'makjennrich: n%d is not a numeric type' , i )
      
    % Is number of samples a valid value?
    elseif  ~ isfinite( n )
      
      error( 'MAK:makjennrich:numeric' , ...
        'makjennrich: n%d is not a valid value i.e. not finite' , i )
      
    % Is number of samples in valid numeric range?
    elseif  n < 2
      
      error( 'MAK:makjennrich:numeric' , ...
        'makjennrich: n%d is less than 2' , i )
      
    end % check number of samples
    
  end % correlation matrices
  
  % Check that size of correlation matrices is the same
  if  p( 1 ) ~= p( 2 )
    
    error( 'MAK:makjennrich:unequal_size' , ...
        'makjennrich: R1 and R2 are not the same size' )
    
  end % R1 R2 same size
  
  % Dimensionality of underlying data
  p = p( 1 ) ;
  
  % Guarantee that all input data is double floating point, for maximum
  % precision
  if  ~ isa( R1 , 'double' ) , R1 = double( R1 ) ; end
  if  ~ isa( n1 , 'double' ) , n1 = double( n1 ) ; end
  if  ~ isa( R2 , 'double' ) , R2 = double( R2 ) ; end
  if  ~ isa( n2 , 'double' ) , n2 = double( n2 ) ; end
  
  
  %%% Jennrich Test %%%
  
  % Degrees of freedom
  df = p .* ( p - 1 ) ./ 2 ;
  
  % Weighted average correlation matrix
  Rbar = ( n1 .* R1  +  n2 .* R2 )  ./  ( n1 + n2 ) ;
  
  % Compute inverse of average correlation matrix once, and use twice
  Rbarinv = pinv( Rbar ) ;
  
  % This is the sum of Kronecker's delta and the element-wise
  % multiplication of Rbar and its inverse
  S = eye( p )  +  Rbar .* Rbarinv ;
  
  % Variable c, used to derive Z
  c = n1 .* n2 ./ ( n1 + n2 ) ;
  
  % Z, which is reminiscent of Linear Discriminant Analysis, is the
  % normalised difference between the correlation matrices
  Z = sqrt( c )  *  Rbarinv  *  ( R1  -  R2 ) ;
  
  % The diagonal of Z
  dgZ = diag( Z ) ;
  
  % At last, we can compute our chi-square value
  chi2 = 0.5 .* trace( Z * Z )  -  dgZ' * pinv( S ) * dgZ ;
  
  % Look up the p-value
  pval = 1 - chi2cdf( chi2 , df ) ;
  
  
end % makjennrich

