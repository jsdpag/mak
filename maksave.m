
function  maksave (  dname  ,  d  ,  vname  ,  varargin  )
% 
% maksave (  dname  ,  d  ,  vname  ,  var1  ,  var2  ,  ...  varN  )
% 
% MET Analysis Toolkit. Save auxiliary data in directory dname for each
% experiment in event struct d, as returned by makload. Cell array of
% strings vname lists the name of each saved variable. Following this,
% there must be an input argument for each name in vname. These arguments
% contain the contents of each variable. Each must be a cell array
% containing the variable data for each experiment in d. Hence, there must
% be numel( vname ) tailing input arguments that each have numel( d )
% cells. The output file for the ith experiment in d( i ) and jth variable
% vname{ j } will contain data var[i]{ j }.
% 
% If there is only one variable, then vname may be string.
% 
% If output directory dname does not exist then it will be created. Output
% file names will have the format:
%
%   <subject_id>.<experiment_id>[.<tag1>.<tag2>.<etc>].mat
% 
% The experiment id is taken to be the minimum found in
% d( i ).experiment_id. If more than one experiment has the same id then an
% error is thrown. Will not write over existing files.
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  
  %%% Constants  %%%
  
  % Make cellfun return cell array , add these input arguments in sequence
  CELFUN = {  'UniformOutput'  ,  false  } ;
  
  
  %%% Check input %%%
  
  % Too few inputs
  if  nargin  <  4
    
    error (  'MAK:maksave:args'  ,  'maksave: too few input arguments'  )
  
  % Output directory
  elseif  ~ isrow (  dname  )  ||  ~ ischar (  dname  )
    
    error (  'MAK:maksave:dname'  ,  'maksave: dname must be a string'  )
    
  % Event data struct
  elseif  ~ isstruct (  d  )
    
    error (  'MAK:maksave:d'  ,  'maksave: d must be a struct'  )
    
  % Variable names
  elseif  ~ iscellstr (  vname  ) 
    
    % vname can be a string if there is one output variable
    if  nargin ~= 4  ||  ~ isrow( vname )  ||  ~ ischar( vname )
      
      if  nargin == 4
    
        error (  'MAK:maksave:vname4'  ,  ...
          'maksave: vname must be a string or cell array of strings'  )
      
      else
        
        error (  'MAK:maksave:vname'  ,  ...
          'maksave: vname must be a cell array of strings'  )
        
      end
      
    end
    
  end % check primary input
  
  % Mismatched number of variable names and output variables
  if  4 < nargin  &&  numel( vname )  ~=  nargin - 3
    
    error (  'MAK:maksave:numvar'  ,  ...
          'maksave: vname must name each output variable'  )
        
  % All output variables must be cell arrays
  elseif  any (  ~ cellfun(  @iscell  ,  varargin  )  )
    
    error (  'MAK:maksave:cellvar'  ,  ...
          'maksave: output variable data must be in cell arrays'  )
        
  % Output variable data for each experiment
  elseif  any(  cellfun(  @( v ) numel( d ) ~= numel( v ) ,  varargin  )  )
    
    error (  'MAK:maksave:expvar'  ,  [ 'maksave: output variable ' , ...
      'data must have an element for each experiment in d' ]  )
        
  end % check output variable arguments
  
  
  %%% Build output file names %%%
  
  % First find experiment id's
  eid = cellfun (  @min  ,  { d.experiment_id }  ) ;
  
  % Must all be unique
  if  numel (  eid  )  ~=  numel ( unique(  eid  ) )
    error (  'MAK:maksave:expid'  ,  [ 'maksave: the minimum ' , ...
      'experiment_id must be unique for each element of d' ]  )
  end
  
  % File name cell array
  fname = cell ( size(  d  ) ) ;
  
  % Experiments
  for  i = 1 : numel (  d  )
    
    % File name , no .mat suffix
    fname{ i } = ...
      sprintf (  '%s.%d'  ,  d( i ).subject_id  ,  eid( i )  ) ;
    
    % There are tags
    if  ~ isempty (  d( i ).header( 1 ).tags  )
      
      % Add them
      fname{ i } = [  fname{ i }  ,  '.'  ,  ...
                      strrep( d( 1 ).header( 1 ).tags , ' , ' , '.' )  ] ;
      
    end
    
    % Full path
    fname{ i } = fullfile (  dname  ,  [ fname{ i } , '.mat' ]  ) ;
    
    % File already exists
    if  exist (  fname{ i }  ,  'file'  )
      error (  'MAK:maksave:exists'  ,  'maksave: file exists %s'  ,  ...
        fname{ i }  )
    end
    
  end % experiments
  
  
  %%% Save data %%%
  
  % Make new directory if it does not exist
  if  ~ exist (  dname  ,  'dir'  )
    
    [ status , message , messageid ] = mkdir (  dname  ) ;
    
    % Failed to create directory
    if  ~ status  ,  error (  messageid  ,  message  )  ,  end
    
  end % new dir
  
  % Experiments
  for  i = 1 : numel (  d  )
    
    % Pair variable names to experiment data
    if  ~ iscellstr (  vname  )
      
      % There is one output variable and vname is a string , known from
      % input checking
      S.( vname ) = varargin{ 1 }{ i } ;
      
    else
      
      % Line up variable names in top row and data in bottom row
      S = [  vname( : )'  ;
             cellfun(  @( v ) v{ i }  ,  varargin  ,  CELFUN{ : }  )  ] ;
      
      % Convert into struct
      S = struct (  S{ : }  ) ;
      
    end % pair var name and data
    
    % Save file
    save (  fname{ i }  ,  '-struct'  ,  'S'  )
    
  end % save data
  
  
end % maksave

