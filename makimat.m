
function  i = makimat (  varargin  )
% 
% i = makimat (  n  )
% i = makimat (  a  )
% i = makimat (  a  ,  b  )
% i = makimat (  ...  ,  '-d'  )
% 
% MET Analysis Kit. Returns the logical index matrix accessing a sub-set of
% the upper-triangular half of a square matrix. Given a scalar integer
% value n > 0, i becomes a n x n logical matrix in which the upper
% triangular portion is set to true, and the rest is set to false. When the
% logical vector a with length n is given then only the upper-triangular
% elements i( row , col ) are true where a( row ) and a( col ) are both
% true. When logical vectors a and b with length n are given then i( row ,
% col ) is true when a( row ) and b( col ) are both true.
% 
% By default, the diagonal of i is false. When the final input argument is
% the string '-d' then the diagonal will be true for input argument n, it
% will contain the vector a if b is not given, or it will contain the
% vector a & b.
% 
% Written by Jackson Smith - May 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Number of input and output arguments
  narginchk  (  1  ,  3  )
  nargoutchk (  0  ,  1  )
  
  % Last argument is a char row vector
  if  ischar (  varargin{ end }  )
    
    % The last input argument is -d flag
    if  strcmp (  varargin{ end }  ,  '-d'  )
    
      % Set diagonal
      dflag = true ;

    % Invalid string
    else
      
      error (  'MAK:makimat:badstring'  ,  ...
        'makimat: input argument %d is invalid string'  ,  nargin  )
      
    end % validate string
    
  % No -d flag
  else
    
    % Do not set diagonal
    dflag = false ;
    
  end % -d flag
  
  % Number of input arguments
  narg = nargin  -  dflag ;
  
  % No input arguments other than -d flag
  if  narg  ==  0
    
    error (  'MAK:makimat:onlydflag'  ,  [ 'makimat: -d flag is the ' , ...
      'only input arg. Must give n, a, or a and b' ]  )
    
  % One input
  elseif  narg  ==  1
    
    % Scalar integer value greater than zero
    if   isscalar (  varargin{ 1 }  )  &&  ...
        isnumeric (  varargin{ 1 }  )  &&  ...
         isfinite (  varargin{ 1 }  )  &&  ...
           isreal (  varargin{ 1 }  )  &&  ...
         0  <  varargin{ 1 }  &&  mod (  varargin{ 1 }  ,  1  )  ==  0
        
       % Retrieve input argument n
       n = varargin{ 1 } ;
       
       % Return empty a and b
       a = [] ;  b = [] ;
       
    % Logical vector
    elseif  islogical (  varargin{ 1 }  )  &&  isvector (  varargin{ 1 }  )
      
      % Retrieve input argument a , point b to a
      a = varargin{ 1 } ;
      b = a ;
      
    % Invalid input
    else
      
      error (  'MAK:makimat:onearg'  ,  [ 'makimat: expecting input ' , ...
        'arg 1 to be scalar integer value > 0 or logical vector' ]  )
      
    end % one input arg
    
  % Two input arguments
  elseif  narg  ==  2
    
    % Two logical vectors
    if  islogical (  varargin{ 1 }  )  &&  ...
          isvector (  varargin{ 1 }  )  &&  ...
            islogical (  varargin{ 2 }  )  &&  isvector (  varargin{ 2 }  )
      
      % Retrieve a and b
      a = varargin{ 1 } ;
      b = varargin{ 2 } ;
      
    % Invalid input
    else
      
      error (  'MAK:makimat:twoarg'  ,  [ 'makimat: expecting input ' , ...
        'args 1 and 2 to be logical vectors' ]  )
      
    end % two input args
    
  end % check input
  
  % Logical vector(s) given as input arguments 
  if  ~ isempty (  a  )
    
    % Need to set n , the size of square matrix
    n = numel (  a  ) ;
    
    % Guarantee that a is column vector ...
    if  ~ iscolumn (  a  )  ,  a = reshape (  a  ,  n  ,  1  ) ;  end
    
    % ... and that b is row vector
    if  ~ isrow (  b  )  ,  b = reshape (  b  ,  1  ,  n  ) ;  end
  
  end % logical vectors given
  
  
  %%% Logical index matrix %%%
  
  % Get the upper-triangular portion set to true
  i = triu (  true( n )  ,  1 - dflag  ) ;
  
  % No argument a given so end now
  if  isempty (  a  )  ,  return  ,  end
  
  % Set upper-triangular elements true where a and b are also true
  i( : ) = i  &  bsxfun (  @and  ,  a  ,  b  ) ;
  
  
end % makimat

