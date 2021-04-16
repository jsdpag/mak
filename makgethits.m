
function  h = makgethits( MCC , H , xy )
% 
% h = makgethits( MCC , H , xy )
% 
% MET Analysis Kit. Tests monocular eye positions in xy against hit regions
% in H. Returns logical true in h for each eye position that is in a hit
% region , false otherwise.
% 
% 
% Input
% 
% MCC - MET controller constants, available from metctrlconst or makmetcon.
% 
% H - Numeric array defines list of hit regions associated with one
%   stimulus link. Each row lists a different hit region ; but eye
%   positions can land in any of these to count as looking at the stimulus.
%   A 6 column array defines a set of circular hit regions, while an 8
%   column array defines rectangular ones.
% 
% xy - Array of monocular eye positions. Each column gives a separate eye
%   position. Row 1 has the x-axis coordinates, and row 2 has y-axis
%   coordinates.
% 
% 
% Output
% 
% h - Logical row vector is true for each eye position that lands in any of
%   the hit regions in H.
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Get hit region constants
  C = MCC.SDEF.ptb.hitregion ;
  
  % How many columns define the hit region?
  ncol = size (  H  ,  2  ) ;
  
  % Number of columns does not match those tested here
  if  all (  ncol  ~=  [ 6 , 8 ]  )
    
    error (  'MAK:makgethits:version'  ,  ...
      'makgethits: version incompatibility , check MET vs MAK versions'  )
    
  % Number of cols not supported by these MET controller constants , raise
  % an error
  elseif  all (  C.ncols  ~=  ncol  )
    
    error (  'MAK:makgethits:ncol'  ,  [ 'makgethits: hit regions ' , ...
      'with %d columns seen but MAK only supports those with:%s\n' , ...
      'Check that MET version matches current version of MAK' ]  , ...
      ncol  ,  sprintf( ' %d' , C.ncols )  )
    
  end % wrong number of columns
  
  % Convert hit region matrix to cell array of row vectors
  H = num2cell (  H  ,  2  ) ;
  
  % Compare eye position to hit region according to hit region shape
  switch  ncol
    
    % Circular hit regions
    case  6  ,  I = cellfun (  @( h ) gethitcirc( C , h , xy )  ,  H  , ...
                      'UniformOutput'  ,  false  ) ;
      
    % Rectangular hit regions
    case  8  ,  I = cellfun (  @( h ) gethitrect( C , h , xy )  ,  H  , ...
                      'UniformOutput'  ,  false  ) ;
    
  end % compare to hit region
  
  % Logical vector for first hit region
  h = I{ 1 } ;
  
  % Combine results of all hit regions into one logical vector using
  % logical OR
  for  i = 2 : numel( I )  ,  h = h  |  I{ i } ; end
  
end % gethits


% Compare monocular eye position against circular hit region
function  i = gethitcirc (  C  ,  h  ,  xy  )
  
  % This hit region is ignored , return false vector
  if  ~ h(  C.sixcol.ignore  )
    i = false (  1  ,  size( xy , 2 )  ) ;
    return
  end
  
  % Eye position relative to centre of circle
  xy( 1 , : ) = xy( 1 , : )  -  h( C.sixcol.xcoord ) ;
  xy( 2 , : ) = xy( 2 , : )  -  h( C.sixcol.ycoord ) ;
  
  % Distance of eye position to centre
  d = sqrt ( sum(  xy .^ 2  ,  1  ) ) ;
  
  % Eye position is within hit region
  i = d  <=  h(  C.sixcol.radius  ) ;
  
end % gethitcirc


% Compare monocular eye position against rectangular hit region
function  i = gethitrect (  C  ,  h  ,  xy  )

  % This hit region is ignored , return false vector
  if  ~ h(  C.eightcol.ignore  )
    i = false (  1  ,  size( xy , 2 )  ) ;
    return
  end
  
  % Eye position relative to centre of rectangle
  xy( 1 , : ) = xy( 1 , : )  -  h( C.eightcol.xcoord ) ;
  xy( 2 , : ) = xy( 2 , : )  -  h( C.eightcol.ycoord ) ;
  
  % Rotation matrix absolute values
  c = cosd ( h(  C.eightcol.rotation  ) ) ;
  s = sind ( h(  C.eightcol.rotation  ) ) ;
  
  % Apply clockwise rotation to eye positions so that the coordinate system
  % is in alignment with the rectangle
  % vector.
  xy = [ c , s ; -s , c ]  *  xy ;
  
  % Take the absolute values to simplify comparison with width and height
  xy = abs (  xy  ) ;
  
  % Check if eye position fell within the hit region
  i = xy( 1 , : )  <=  h( C.eightcol.width  ) / 2  &  ...
      xy( 2 , : )  <=  h( C.eightcol.height ) / 2 ;
  
end % gethitrect

