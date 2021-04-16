
function  J = makconnstrength (  E  ,  n  )
% 
% J = makconnstrength (  E  ,  n  )
% 
% MET Analysis Kit, pre-processing. Returns connection-strength matrix J,
% computed from raw interface-energy matrix E and per-cluster spike count
% vector n.
% 
% References:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% UltraMegaSort2000, https://neurophysics.ucsd.edu/software.php
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  % Make sure that spike-count vector is double floating point
  if  ~ isa (  n  ,  'double'  )  ,  n = double (  n  ) ;  end
  
  % Number of clusters
  cnum = numel (  n  ) ;
  
  % Linear index of diagonal elements
  i = find ( eye(  cnum  ) ) ;
  
  % Compute normalised interface energy
  En = n'  *  n ;
  En( i ) = ( En( i )  -  n' )  /  2 ;
  En = E  ./  En ;

  % Compute connection strengths
  J = 2  *  En  ./  ...
    (  repmat( En( i ) , 1 , cnum )  +  repmat( En( i )' , cnum , 1 )  ) ;
  
  % When numerical limitations of the computer cause within-cluster energy
  % to be zero then we will get 0 / 0 = NaN. A sensible thing seems to set
  % the self energy for such clusters to 1.
  J(  i( E( i ) == 0 )  ) = 1 ;
  
end % makconnstrength

