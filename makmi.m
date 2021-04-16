
function  [ I , Is , pval ] = makmi (  X  ,  Y  ,  rep  ,  varargin  )
% 
% [ I , Is , pval ] = makmi (  X  ,  Y  ,  rep  )
%        ... = makmi (  ...  ,  nbins  )
%        ... = makmi (  ...  ,  Xedges  ,  Yedges  )
%        ... = makmi (  ...  ,  Name  ,  Value  )
% 
% MET Analysis Kit. Computes the empirical mutual information between input
% X and output Y, both containing numeric or logical arrays. X is accessed
% by linear indexing as though it were a vector. If X has N elements then Y
% must be either an N-element vector or an array with size N along
% dimension 1 i.e. multi-dimensional Y must have as many rows as there are
% elements in X. In either case, the elements ( Y is vector ) or rows ( Y
% is array ) of Y must be in register with the elements of X. If Y is
% multi-dimensional then makmi computes mutual information separately for
% each column of Y.
% 
% Be aware that X or Y can be integer identifiers that categorise the
% original values in some way. For example, each element of Y might be an
% integer that uniquely identifies one of the 2 ^ T possible binary spike
% trains over T time bins that was observed on that trial.
% 
% I contains the mutual information between X and Y that is corrected for
% bias by subtracting an estimated shuffle-correcting value, returned in
% Is. Hence, the raw and biased mutual information is I + Is. The shuffle
% correction is done by randomly permuting X and recomputing the mutual
% information with Y. This is repeated rep times ( default 30 ) and the
% average of the resulting distribution is used for Is ; rep is optional ,
% or the default value can be used by providing an empty matrix [ ] ;
% shuffling can be skipped altogether by setting rep to zero. If Y is a
% vector then both I and Is are scalars. If Y is multi-dimensional with
% size s then both I and Is will have size s( 2 : end ), where the ith
% element of I or Is corresponds to the ith column of Y i.e. Y( : , i ),
% where i is a linear index. pval is the estimated significance of the
% empirical value in I. It is computed from the distribution of shuffled
% information values using a one-tailed test such that pval is the fraction
% of shuffled values that exceed or equal I + Is. pval has the same size as
% I and Is. If rep is zero then Is will contain zeros and pval will contain
% ones.
% 
% Because the empirical joint probability distribution of X and Y is
% obtained using histcounts2, the way that the data is binned can be
% controlled using additional optional inputs nbins or *edges, and Name /
% Value pairs. See doc histcounts2 for details. Note that 'Normalization'
% is silently ignored to guarantee that histcounts2 always returns the
% count.
% 
% For example, X might be the orientation of a sinusoidal grating on each
% of N trials and Y might be the firing rate of a visual neurone. On the
% ith trial, the grating had orientation X( i ) and the neurone fired at a
% rate of Y( i ) in response. Using the default binning algorithm of
% histcounts2, the mutual information is I = makmi (  X  ,  Y  ).
% 
% Note that makmi uses parfor loops to compute shuffle-correctors.
% 
% 
% References:
% 
%   Hatsopoulos NG, et al. 1998. Information about movement direction
%     obtained from synchronous activity of motor cortical neurons. Proc
%     Natl Acad Sci U S A 95(26): 15706-15711.
%
%   Rieke F, Warland D, de Ruyter van Steveninck R, Bialek W. 1997. Spikes:
%     Exploring the neural code. The MIT Press, London England. ISBN
%     0262181746.
%   
%   Wallisch, Pascal, et al. MATLAB for Neuroscientists : An Introduction
%     to Scientific Computing in MATLAB, Elsevier Science & Technology,
%     2014. ISBN 9780123838377.
% 
% 
% Written by Jackson Smith - May 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input %%%
  
  % Is X numeric?
  if  ~ isnumeric (  X  )  &&  ~ islogical (  X  )
    
    error (  'MAK:makmi:X'  ,  'makmi: X must be numeric or logical'  )
    
  end % X numeric
  
  % Get the number of elements in X
  N = numel (  X  ) ;
  
  % Is X empty?
  if  N  ==  0
    
    % Return empty matrices
    I = [] ;  Is = [] ;
    return
    
  % X is not a column vector , make it so
  elseif  ~ iscolumn (  X  )
    
    X = reshape (  X  ,  N  ,  1  ) ;
    
  end % X is empty or not column vect
  
  % Is Y numeric?
  if  ~ isnumeric (  Y  )  &&  ~ islogical (  Y  )
    
    error (  'MAK:makmi:Ynumeric'  ,  ...
      'makmi: Y must be numeric or logical'  )
    
  end % Y numeric
  
  % The size of Y in each dimension
  sy = size (  Y  ) ;
  
  % Y is a vector
  if  isvector (  Y  )
    
    % Does it have as many elements as X?
    if  numel (  Y  )  ~=  N
      
      error (  'MAK:makmi:Ynumel'  ,  ...
        'makmi: vector Y must have as many elements as X'  )
      
    % Y is not a column vector , make it so
    elseif  ~ iscolumn (  Y  )
      
      Y = reshape (  Y  ,  N  ,  1  ) ;
      sy = [ N , 1 ] ;
      
    end % numel Y matches X
    
  % Y is not a vector
  else
    
    % Does Y have as many rows as there are elements in X?
    if  sy( 1 )  ~=  N
      
      error (  'MAK:makmi:Ynumrows'  ,  [ 'makmi: array Y must have ' , ...
        'as many rows as there are elements in X' ]  )
      
    end % size Y dim 1 matches numel X
    
  end % size of Y
  
  % Columns of Y
  Ncol = prod (  sy( 2 : end )  ) ;
  
  % rep is missing , use default value
  if  nargin  <  3  ||  isempty (  rep  )  ,  rep = 30 ;  end
  
  % Rep is a scalar integer value greater than or equal to zero
  if  ~ isscalar (  rep  )  ||  ~ isnumeric (  rep  )  ||  ...
        ~ isreal (  rep  )  ||  ~ isfinite (  rep  )  ||  ...
          rep  <  0  ||  mod (  rep  ,  1  )
    
    error (  'MAK:makmi:rep'  ,  ...
      'makmi: rep must be a scalar real integer value of zero or more '  )
    
  end % check rep
  
  % Look for histcounts2 'Normalization' Name/Value pairs
  i = find ( strcmp(  varargin  ,  'Normalization'  ) ) ;
  
  % And get rid of them , we must guarantee that the count is returned
  varargin( [ i , i + 1 ] ) = [] ;
  
  
  %%% Mutual information %%%
  
  % Allocate output arrays
    I  = zeros (  1  ,  Ncol  ) ;
    Is = zeros (  1  ,  Ncol  ) ;
  pval =  ones (  1  ,  Ncol  ) ;
  
  % Cell array for holding temporary bin edges
  edges = cell (  1  ,  2  ) ;
  
  % Columns of Y
  for  i = 1 : Ncol
    
    % Grab column of Y
    Yc = Y( : , i ) ;
    
    % Determine what the bin edges are JUST ONCE, then pass these into the
    % parfor loop. We also take the opportunity to calculate empirical Pxy.
    [ ePxy , edges{ : } ] = histcounts2 (  X  ,  Yc  ,  varargin{ : }  ) ;
      
      % Must normalise after getting output
      ePxy = ePxy  /  N ;
    
    % Calculate the raw empirical mutual information
    I( i ) = mi (  X  ,  Yc  ,  N  ,  false  ,  false  ,  { ePxy }  ) ;
    
    % No shuffling , go to next column
    if  ~ rep  ,  continue  ,  end
    
    % Build distribution of mutual information from randomly permuted data
    parfor  j = 1 : rep
      
      smi( j ) = mi (  X  ,  Yc  ,  N  ,  true  ,  true  ,  edges  ) ;
      
    end % shuffled mutual info
    
    % Take the average shuffled mutual information as the correction term
    Is( i ) = mean (  smi  ) ;
    
    % Count the number of shuffled information values greater than or equal
    % to the empirical mutual information
    pval( i ) = sum (  I( i )  <=  smi  ) ;
    
  end % Y cols
  
  % Shuffle correction was done
  if  rep
    
    % Subtract estimate of bias from mutual information
    I( : ) = I  -  Is ;
    
    % Estimate significance of empirical mutual information , divide counts
    % by total number
    pval( : ) = pval  ./  rep ;
  
  end % apply shuffle correction
  
  % Reshape output to match upper dimensions of input Y
  I    = reshape (  I     ,  [ sy( 2 : end ) , 1 ]  ) ;
  Is   = reshape (  Is    ,  [ sy( 2 : end ) , 1 ]  ) ;
  pval = reshape (  pval  ,  [ sy( 2 : end ) , 1 ]  ) ;
  
  
end % makmi


%%% Subroutine %%%

% Computes empirical mutual information between vectors x and y , each with
% N elements. rp is true if x should first be randomly permuted. Binning
% instructions are provided in variable input argument cell array vin if
% dohc is true, otherwise vin contains Pxy.
function  i = mi (  x  ,  y  ,  N  ,  rp  ,  dohc  ,  vin  )
  
  % Randomly permute x
  if  rp  ,  x = x(  randperm(  N  )  ) ;  end
  
  % Do we need to do 2D binning?
  if  dohc
    
    % Two-dimensional binning divided by the number of elements gives the
    % empirical joint probability distribution of x and y
    Pxy = histcounts2 (  x  ,  y  ,  vin{ : }  )  /  N ;
  
  % Nope
  else
    
    % Just get Pxy from input argument
    Pxy = vin{ 1 } ;
    
  end % To bin or not to bin?
  
  % From this we get the marginal distributions
  Px = sum (  Pxy  ,  2  ) ;
  Py = sum (  Pxy  ,  1  ) ;
  
  % The joint probability is normalised by the product of the marginal
  % probabilities before taking the logarithm. Compute the product of
  % marginal probabilities for all combinations of x and y.
  PxPy = bsxfun (  @times  ,  Px  ,  Py  ) ;
  
  % Reshape Pxy and PxPy into column vectors , for convenience. Note that
  % joint probabilities and the product of marginal probabilities are still
  % in register with respect to each data bin.
  Pxy = Pxy( : ) ;  PxPy = PxPy( : ) ;
  
  % Finally , compute the mutual information
  i = sum (  Pxy  .*  log2 (  Pxy  ./  PxPy  )  ,  'omitnan'  ) ;
  
end % mi

