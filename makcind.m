
function  [ i1 , i2 , j ] = makcind (  c1  ,  c2  ,  cnum  )
% 
% [ i1 , i2 , j ] = makcind (  c1  ,  c2  ,  cnum  )
% 
% MET Analysis Kit, pre-processing. Returns the linear indices for a
% cnum x cnum interface energy matrix that are needed to re-compute the
% interface energy following a merger of cluster c1 and cluster c2, where
% c1 < c2. i1 has the indices for all energies between c1 and every other
% cluster, except c2. i2 has the indices for all energies between c2 and
% every other cluster, except c1. j is a sorted list of all cluster indices
% except for c1 and c2. i1 and i2 both define 'L' shaped regions in the
% interface energy matrix.
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Determine the subscript row and column indices needed to find i1
  r = [  1 : c1 - 1  ,  repmat( c1 , 1 , cnum - c1 - 1 )  ] ;
  c = [  repmat( c1 , 1 , c1 - 1 )  ,  ...
       [ c1 + 1 : c2 - 1 , c2 + 1 : cnum ]  ] ;
  
  % Get equivalent linear indices
  i1 = sub2ind (  [ cnum , cnum ]  ,  r  ,  c  ) ;
  
  % Do this again for i2
  r = [  [ 1 : c1 - 1 , c1 + 1 : c2 - 1 ]  ,  ...
    repmat( c2 , 1 , cnum - c2 )  ] ;
  c = [  repmat( c2 , 1 , c2 - 2 )  ,  c2 + 1 : cnum  ] ;
  i2 = sub2ind (  [ cnum , cnum ]  ,  r  ,  c  ) ;
  
  % List of cluster indices lacking c1 and c2
  j = [  1 : c1 - 1  ,  c1 + 1 : c2 - 1  ,  c2 + 1 : cnum  ] ;
  
end % makcind

