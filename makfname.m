
function  s = makfname (  d  ,  short  )
% 
% s = makfname (  d  ,  short  )
% 
% MET Analysis Toolkit. Returns a file name in s ( without path ) from
% experiment data struct d , as returned by makload. The format is:
% 
%   < subject_id >.< experiment_id 1 >[.< eid 2 >.< ... >]
%     [.< tag 1 >.< tag 2 >.< ... >].mat
% 
% If d consists of multiple elements then the unique set of experiment id's
% are chained together. Subject id and tags are always obtained from
% d( 1 ). When short is non-zero and there are three or more experiments
% ( elements ) in d then the format compresses to:
% 
%   < subject_id >.< eid 1 >.to.< eid N >...
% 
% If short is not provided then the first format is used.
% 
% This is suitable for naming files that store results grouped over
% multiple experiments.
% 
% Written by Jackson Smith - DPAG
% 
  
  % Check input and output argument number
  narginchk  (  1  ,  2  )
  nargoutchk (  0  ,  1  )
  
  % Input must be a struct
  if  ~ isstruct (  d  )
    
    error (  'MAK:makfname:inputd'  ,  ...
      'makfname: input d must be a struct array'  )
    
  % short is provided
  elseif  2  <=  nargin
    
    if  ~ isscalar (  short  )  ||  ...
        ~ (  isnumeric(  short  )  ||  islogical(  short  )  )
    
      error (  'MAK:makfname:inputshort'  ,  ...
        'makfname: input short must be a scalar numeric or logical'  )
    
    end
    
  % short not provided , take default
  else
    
    short = 0 ;
    
  end % input check
  
  % Subject id
  subid = d( 1 ).subject_id ;
  
  % Tag string
  if  ~ isempty (  d( 1 ).header( 1 ).tags  )
    tags = [  '.'  ,  strrep( d( 1 ).header( 1 ).tags , ' , ' , '.' )  ] ;
  else
    tags = '' ;
  end
  
  % Unique set of experiment ids
  eid = unique (  [ d.experiment_id ]  ) ;
  
  % Use short format eid chain
  if  2  <  numel (  eid  )  &&  short
    
    eid = sprintf (  '.%d.to.%d'  ,  eid( [ 1 , end ] )  ) ;
    
  % Long format eid chain
  else
    
    eid = sprintf (  '.%d'  ,  eid  ) ;
  
  end % eid chain
  
  % File name string
  s = [  subid  ,  eid  ,  tags  ,  '.mat'  ] ;
  
end % makfname

