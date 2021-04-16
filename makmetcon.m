
function  [ MC , MCC ] = makmetcon %#ok
% 
% [ MC , MCC ] = makmetcon
% 
% MET Analysis Kit. Attempts to load and return MET compile time constants
% and MET controller constants that are packaged with MAK.
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Get directory containing MAK files
  dname = fileparts ( which(  'makmetcon'  ) ) ;
  
  % Constant .mat file names
  fmc  = fullfile (  dname  ,  'MC.mat'   ) ;
  fmcc = fullfile (  dname  ,  'MCC.mat'  ) ;
  
  % Make sure that files exist
  if  ~ exist (  fmc  ,  'file'  )
    
    error (  'MAK:makmetcon:MC'  ,  'makmetcon:can''t find MC.mat'  )
    
  elseif  ~ exist (  fmcc  ,  'file'  )
    
    error (  'MAK:makmetcon:MCC'  ,  'makmetcon:can''t find MCC.mat'  )
    
  end % check existence of files
  
  % Load constants
  load (  fmc  ,  'MC'  )
  load (  fmcc ,  'MCC' )
  
end % makmetcon

