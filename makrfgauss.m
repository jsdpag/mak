
function  [ C , E , R2 , coord , SEL , MAP ] = makrfgauss( RFMAP , VARNAM )
% 
% [ C , E , R2 , coord , SEL , MAP ] = makrfgauss( RFMAP , VARNAM )
% 
% MET Analysis Kit.
% 
% Simple GUI tool for user-guided receptive field (RF) quantification. The
% user provides manual guidance in fitting a 2D isotropic Gaussian to RF
% mapping data for each electrode. First, the user chooses which RF mapping
% data set(s) to include, then provides starting baseline, centre and width
% parameters. Least-squares linear regression fits are then generated.
% 
% RFMAP is a string naming the makprep output file containing RF mapping
% data. It is assumed that no spike sorting was done and that the d.spike
% field is as specified in doc makprep. All spikes from the analysis epoch
% are counted; it is recommended that pre-processing of RF mapping sessions
% uses an epoch that covers only the response to the mapping stimulus. Only
% the multiunit spiking will be used to determine RF positions, under the
% assumption that units on the same electrode will have very similar
% retinal topography.
% 
% If more than one RF mapping session is available then RFMAP can be a cell
% array of strings naming the pre-processed files. There will be an
% additional step in the manual guidance.
% 
% For each electrode, the user will be asked to do the following:
% 
%   1) Select RF mapping sessions. (skipped if a single session is given)
%   2) Specify baseline values. (For each RF mapping session)
%   3) Choose centre. (One centre applies to all RF mapping sessions)
%   4) Set width. (Applies to all)
%   5) Verify least-square's fit.
% 
% Basic controls of GUI. In step 1), an axes is presented for each RF
% mapping session. RF mapping average firing rates are shown in greyscale
% (white to black is zero to max rate). Left-click each axes object to
% include that RF mapping session. Left click again to de-select. Shift-
% click to select all; shift-click again to de-select all. Selection for
% the previous electrode is set by default for the next. In step 2), left
% click points on each selected RF map to specify baseline values. Data at
% selected baseline coordinates are not used for the least-squares fitting.
% Either right-click or control-click an axis to initiate sampling of two
% coordinates; all RF map values within range of these samples are
% selected. Step 3), left click to set the centre point. Step 4), left
% click to set the radius of a circle around the centre point. In step 5),
% a proofing figure is made, showing the best-fitting Gaussians for visual
% comparison against the empirical RF maps; the maps are updated to show
% the best fitting centre and width. Hit <Enter> or close the proofing
% figure to accept the fit. Hit <Esc> to reject the fit and try again with
% the same electrode from step 1). In steps 1 to 4, hit <Enter> to advance
% to the next step, or <Esc> to clear current selection.
% 
% If the main figure is ever closed then the program aborts and returns
% only empty matrices in all outputs. One may roll back to an earlier
% electrode and then re-fit the following electrodes with by typing
% Shift-Escape and then giving the index of the desired electrode.
% 
% The median baseline is then subtracted from each session, and the
% resulting peak is found within range of the given centre and width. Each
% session is normalised by that peak so that all values vary in the range
% -1 to +1. This step is used to account for differences in noise level and
% signal quality between recording days.
% 
% A single set of coefficients are found that optimise the least-squares
% fit of a 2D isotropic Gaussian to all of the normalised RF mapping data.
% 
% Optional input VARNAM can be provided to name the pair of task variables
% that contain the x- and y-axis coordinate of the mapping stimulus. This
% will be a two element cell array of strings. If not provided then VARNAM
% is assumed to be { 'TargXCoord' , 'TargYCoord' }. We make the terrible
% assumption that the SAME variable names are used in all RF mapping
% sessions named in RFMAP. It is required that the x-coordinate variable is
% in VARNAM{ 1 } and that the y-coordinate variable is in VARNAM{ 2 }.
% 
% Returns output arguments:
% 
%   C - Electrode x Coefficient matrix. Rows are indexed by electrode.
%     Columns contain Guassian coefficients, ordered [ x , y , w , a , b ].
%     The coordinate (x,y) is the centre of the RF, w is one standard
%     deviation, a is the amplitude, and b is the baseline. x, y, and w are
%     in degrees of visual field. a is in arbitrary units, being fit to
%     normalised data. The coefficients can be used to recreate the shape
%     of the RF at point (a,b) with:
%     
%       f(a,b) = a * exp( -0.5 * z * C * z' )  +  b
% 
%         where z = [ a , b ] - [ x , y ] and
%               C = 1 / w ^ 2 * [ 1 , 0 ; 0 , 1 ]
%   
%   E - Column vector with electrode ID values, in register with rows of C.
%   
%   R2 - Column vector in register with rows of C and elements of E.
%     Contains the coefficient of determination i.e. R^2 for each gaussian
%     fit.
%   
%   coord - struct vector containing the marginal x- and y-axis coordinates
%     of each RF map given in RFMAP. For use with SEL and MAP.
%   
%   SEL - Electrode x RF map logical matrix, each row contains a logical
%     index vector flagging which RF maps were used for each electrode.
%   
%   MAP - Electrode x RF maps cell array, each cell contains a matrix of
%     the average firing rate in response to each mapping stimulus
%     location. Rows index x-axis coordinates, columns index over y-axis
%     coordinates, ascending in each case. All in register with coord.
% 
% Written by Jackson Smith - April 2020 - DPAG, University of Oxford
% 
  
  
  %%% Constants %%%
  
  % Number of gaussian coefficients to fit
  n.coef = 5 ;
  
  % Max number of electrodes
  n.maxein = 128 ;
  
  % Figure properties
  CONST.fig = { 'MenuBar' , 'none' , 'DockControls' , 'off' , ...
    'Colormap' , flip( gray( 256 ) , 1 ) , 'Tag' , 'mainfig' } ;
  
  % Axes properties
  CONST.ax = { 'TickDir' , 'out' , 'FontSize' , 12 , 'LineWidth' , 1 , ...
    'NextPlot' , 'add' } ;
  
  % Scatter starting properties
  CONST.sc = { [] , [] , 32 , 'filled' , 'MarkerEdgeColor' , '#77AC30' ,...
    'MarkerFaceColor' , '#77AC30' , 'PickableParts' , 'none' , ...
      'Tag' , 'baseline' } ;
  
  % Line starting properties
  CONST.ln = { NaN , NaN , 'Marker' , '+' , 'LineWidth' , 2 , ...
    'LineStyle' , 'none' , 'Color' , '#EDB120' , 'Tag' , 'centre' , ...
      'MarkerSize' , 8 , 'PickableParts' , 'none' } ;
  
  % Rectangle starting properties
  CONST.rt = { 'Position' , zeros( 1 , 4 ) , 'LineWidth' , 2 , ...
    'EdgeColor' , '#EDB120' , 'Curvature' , [ 1 , 1 ] , ...
      'PickableParts' , 'none' , 'Tag' , 'width' } ;
  
  % Point to root graphics object
  CONST.rgo = groot ;
  
  % lsqcurvefit options, turn off messages
  CONST.lop = optimoptions( 'lsqcurvefit' , 'Display' , 'off' ) ;
  
  
  %%% Check input %%%
  
  % Make sure number of input and output arguments fits within valid range
  narginchk ( 1 , 2 )
  nargoutchk( 0 , 6 )
  
  % VARNAM not provided
  if  nargin < 2
    
    % Default
    VARNAM = { 'TargXCoord' , 'TargYCoord' } ;
    
  % VARNAM is provided, make sure it is a cell array of two strings
  elseif  ~ iscellstr( VARNAM )  ||  numel( VARNAM ) ~= 2  %#ok
    
    error( 'MAK:makrfgauss:VARNAM' , [ 'makrfgauss: VARNAM ' , ...
        'must be a 1 x 2 cell array of strings' ] )
    
  end % VARNAM

  % RFMAP is a string
  if  ischar( RFMAP )

    % Then pack it into a cell array. Guarantee row vector.
    RFMAP = { RFMAP( : )' } ;

  end % RFMAP is string

  % Is RFMAP a cell array of strings?
  if  ~ iscellstr( RFMAP )

    error( 'MAK:makrfgauss:RFMAP' , ...
      'makrfgauss: RFMAP must be a string or a cell array of strings' )

  end % cell array of strings?
  
  % Number of file names
  n.rfmap = numel( RFMAP ) ;
  
  % Allocate cell array for loading d structs
  d = cell( 1 , n.rfmap ) ;
  
  % Also, alocate mapping coordinates
  coord.x = cell( 1 , n.rfmap ) ;
  coord.y = cell( 1 , n.rfmap ) ;

  % Check that each string names a file containing a variable called 'd'
  for  i = 1 : n.rfmap
    
    % Does this file exist?
    if  ~ exist( RFMAP{ i } , 'file' )
      
      error( 'MAK:makrfgauss:RFMAP_exist' , [ 'makrfgauss: file name ' ,...
        '%d in RFMAP is not a valid file name' ] , i )
      
    end % file exists?
    
    % The file exists, now look for variable called 'd'
    S = whos( '-file' , RFMAP{ i } , 'd' ) ;
    
    % Struct is empty if variable not found
    if  isempty( S )
      
      error( 'MAK:makrfgauss:RFMAP_nod' , [ 'makrfgauss: file ' , ...
        '%d in RFMAP does not contain variable ''d''' ] , i )
      
    end % 'd' not found
    
    % 'd' was found, now load it
    d{ i } = load( RFMAP{ i } , 'd' ) ;
    d{ i } = d{ i }.d ;
    
    % Look for variable names in VARNAM
    if  ~ all( isfield( d{ i }.sd( 1 ).var , VARNAM ) )
      
      error( 'MAK:makrfgauss:VARNAM_missing' , [ 'makrfgauss: file ' , ...
        '%d in RFMAP does not contain both variables %s and %s' ] , i , ...
          VARNAM{ : } )
      
    end % variables not found
    
    % Get unique set of mapping coordinates
    coord.x{ i } = unique( d{ i }.sd( 1 ).var.( VARNAM{ 1 } ).value ) ;
    coord.y{ i } = unique( d{ i }.sd( 1 ).var.( VARNAM{ 2 } ).value ) ;
    
  end % file names
  
  % Collapse coord and d into struct row vectors
  coord = struct( 'x' , coord.x , 'y' , coord.y ) ;
  d = [ d{ : } ] ;
  
  % Count number of marginal coordinates
  n.coord.x = arrayfun( @( c ) numel( c.x ) , coord ) ;
  n.coord.y = arrayfun( @( c ) numel( c.y ) , coord ) ;
  
  
  %%% Preparation %%%
  
  % List of all electrodes that had at least 1 spike
  E = unique( [ d.electrodes ] )' ;
  
  % Number of electrodes
  n.eid = numel( E ) ;
  
  % Allocate remaining output
    C = zeros( n.eid , n.coef  ) ;
   R2 = zeros( n.eid ,       1 ) ;
  SEL = false( n.eid , n.rfmap ) ;
  MAP =  cell( n.eid , n.rfmap ) ;
  
  % Determine number of rows and columns in presenting RF maps
  if  n.rfmap == 1
    n.ax.col = 1 ;
    n.ax.row = 1 ;
  else
    n.ax.col = 2 ;
    n.ax.row = ceil( n.rfmap ./ 2 ) ;
  end
  
  % Create figure object
  h.f = figure( CONST.fig{ : } , 'KeyPressFcn' , @fkeypress ) ;
  
  % If number of RF mapping sessions is more than 4, then stretch figure
  % from top to bottom of screen
  if  4  <  n.rfmap
    
    h.f.OuterPosition( [ 2 , 4 ] ) = CONST.rgo.ScreenSize( [ 2 , 4 ] ) ;
    
  end
  
  % Initialise logical vector indicating which RF maps are selected. Start
  % with none selected. Vector in register with RFMAP.
  h.f.UserData.sel = false( n.rfmap , 1 ) ;
  
  % Valid RF maps have data for current electrode
  h.f.UserData.val = false( n.rfmap , 1 ) ;
  
  % Initialise figure keypress name
  h.f.UserData.key = '' ;
  
  % Maximum number of electrodes
  h.f.UserData.maxein = n.maxein ;
  
  % Manual reset of eind. If non-zero, it is the new eind. If zero then
  % flag is low, no new eind given.
  h.f.UserData.new_eind = 0 ;
  
  
  %%% Least-squares fitting %%%
  
  % Electrode index, runs from 1 to n.eid
  eind = 1 ;
  
  % Electrode loop
  while  eind <= n.eid
    
    % Calculate RF map for this electrode, in each session.
    M = calcrfmap( VARNAM , coord , d , E( eind ) ) ;
    
      % Copy into output
      MAP( eind , : ) = M ;
    
    % Check which maps are valid
    h.f.UserData.val( : ) = ~ cellfun( @isempty , M ) ;
    
    % This electrode is not represented in any mapping session
    if  ~ any( h.f.UserData.val )
      
      % Flag as invalid electrode
      C( eind , : ) = NaN ;
      
      % Next trode
      eind = eind  +  1 ;
      continue
      
    end % no mapping data for electrode
    
    % Make sure to de-select RF maps that lack any data
    h.f.UserData.sel( : ) = h.f.UserData.sel  &  h.f.UserData.val ;
    
    
    %- STEP 1) Select RF maps -%
    
    % Remember which step this is
    if  fignam( h.f , 1 , E , eind , n.eid , 'Chose RF maps' )
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    end

    % Initialise maps
    for  i = 1 : n.rfmap

      % Create axes
      h.ax = subplot( n.ax.row , n.ax.col , i , 'UserData' , i , ...
        CONST.ax{ : } , 'YDir' , 'Normal' , 'ButtonDownFcn' , @a1click );

      axis  equal tight
      colorbar
      
      % Turn off the AxesToolbar object
      h.ax.Toolbar.Visible = 'off' ;

      % Show map, if it exists
      if  ~ isempty( M{ i } )
        imagesc( h.ax , coord( i ).x , coord( i ).y , M{ i }' , ...
          'PickableParts' , 'none' )
        axis( [ xlim( h.ax ) , ylim( h.ax ) ] )
      end

      % Indicate selection status
      if  h.f.UserData.sel( i )  &&  h.f.UserData.val( i )
        title( h.ax , 'SELECTED' )
      end

    end % init maps
    
    % Step 1) select RF maps, but only if there is more than one map
    if  1  <  n.rfmap
      
      % Wait for user to hit <Enter> following manual selections
      uiwait( h.f )
      
      % New electrode index given
      if  ishandle( h.f ) && h.f.UserData.new_eind
        eind = h.f.UserData.new_eind ; % Copy to master eind variable
        h.f.UserData.new_eind = 0 ;    % Lower flag
        delete( h.f.Children )         % Kill axes
        continue                       % New iteration of electrode loop
      end
      
    else
      
      % Make sure that the single map is selected
      h.f.UserData.sel( : ) = 1 ;
      title( h.ax , 'SELECTED' )
      
    end % Step 1)
    
    % Get rid of Axes mouse click function
    h.ax = findobj( h.f , 'Type' , 'axes' ) ;
    set( h.ax , 'ButtonDownFcn' , [ ] )
    
    
    %- STEP 2) Choose baseline values -%
    
    if  fignam( h.f , 2 , E , eind , n.eid , 'Baseline' )
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    end
    
    % Copy selection into output
    SEL( eind , : ) = h.f.UserData.sel ;
    
    % RF maps
    for  i = 1 : n.rfmap
      
      % Not selected, continue to next
      if  ~ h.f.UserData.sel( i )  ,  continue  ,  end
      
      % Locate axes with this map index
      h.ax = findobj( h.f , 'UserData' , i ) ;
      
      % Create scatter object to show selected points
      scatter( h.ax , CONST.sc{ : } )
      
      % Update mouse-click callback
      h.ax.ButtonDownFcn = @( h , ~ ) a2click( h , coord( i ) , M{ i } ) ;
      
    end % RF maps
    
    % Wait for user to select points and hit <Enter>
    uiwait( h.f )
    
    % New electrode index given
    if  ishandle( h.f ) && h.f.UserData.new_eind
      eind = h.f.UserData.new_eind ; % Copy to master eind variable
      h.f.UserData.new_eind = 0 ;    % Lower flag
      delete( h.f.Children )         % Kill axes
      continue                       % New iteration of electrode loop
    end
    
    % Get rid of Axes mouse click function
    h.ax = findobj( h.f , 'Type' , 'axes' ) ;
    set( h.ax , 'ButtonDownFcn' , [ ] )
    
    
    %- STEP 3) Choose centre -%
    
    if  fignam( h.f , 3 , E , eind , n.eid , 'Centre' )
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    end
    
    % RF maps
    for  i = 1 : n.rfmap
      
      % Not selected, continue to next
      if  ~ h.f.UserData.sel( i )  ,  continue  ,  end
      
      % Locate axes with this map index
      h.ax = findobj( h.f , 'UserData' , i ) ;
      
      % Create line object to show centre of Gaussian
      h.c = line( h.ax , CONST.ln{ : } ) ;
      
        % Now empty the thing
        set( h.c , 'XData' , [] , 'YData' , [] )
      
      % Update mouse-click callback
      h.ax.ButtonDownFcn = @a3click ;
      
    end % RF maps
    
    % Wait for user to select centre and hit <Enter>
    uiwait( h.f )
    
    % New electrode index given
    if  ishandle( h.f ) && h.f.UserData.new_eind
      eind = h.f.UserData.new_eind ; % Copy to master eind variable
      h.f.UserData.new_eind = 0 ;    % Lower flag
      delete( h.f.Children )         % Kill axes
      continue                       % New iteration of electrode loop
    end
    
    % Get rid of Axes mouse click function
    h.ax = findobj( h.f , 'Type' , 'axes' ) ;
    set( h.ax , 'ButtonDownFcn' , [ ] )
    
    
    %- STEP 4) Choose centre -%
    
    if  fignam( h.f , 4 , E , eind , n.eid , 'Width' )
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    end
    
    % RF maps
    for  i = 1 : n.rfmap
      
      % Not selected, continue to next
      if  ~ h.f.UserData.sel( i )  ,  continue  ,  end
      
      % Locate axes with this map index
      h.ax = findobj( h.f , 'UserData' , i ) ;
      
      % Create rectangle object to show 1 standard dev of Gaussian
      rectangle( h.ax , CONST.rt{ : } )
      
      % Update mouse-click callback
      h.ax.ButtonDownFcn = @a4click ;
      
    end % RF maps
    
    % Wait for user to select centre and hit <Enter>
    uiwait( h.f )
    
    % New electrode index given
    if  ishandle( h.f ) && h.f.UserData.new_eind
      eind = h.f.UserData.new_eind ; % Copy to master eind variable
      h.f.UserData.new_eind = 0 ;    % Lower flag
      delete( h.f.Children )         % Kill axes
      continue                       % New iteration of electrode loop
    end
    
    % Get rid of Axes mouse click function
    h.ax = findobj( h.f , 'Type' , 'axes' ) ;
    set( h.ax , 'ButtonDownFcn' , [ ] )
    
    
    %- STEP 5) Validate gaussian fit -%
    
    if  fignam( h.f , 5 , E , eind , n.eid , 'Validate' )
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    end
    
    % Initialise approval flag. If down then the user has not accepted this
    % fit and it will be done again. If up then the user has accepted the
    % fit. We raise the flag in case the user closes the proofing figure,
    % to indicate approval of the fit.
    h.f.UserData.appflg = true ;
    
    % Find centre and width
    h.c = findobj( h.f , 'Tag' , 'centre' ) ;
    h.w = findobj( h.f , 'Tag' ,  'width' ) ;
    
      % They all have same values, get first of each set
      h.c( 2 : end ) = [ ] ;
      h.w( 2 : end ) = [ ] ;
    
    % Centre x- and y-coordinate
    x = h.c.XData ;  y = h.c.YData ;
    
    % Width
    w = x - h.w.Position( 1 ) ;
    
    % Allocate X and Y coordinate accumulators
    X = cell( size( M ) ) ;
    Y = cell( size( M ) ) ;
    
    % Find normalisation value for each RF map, nan flags cases where
    % selected centre and width fall off the RF map
    mval = zeros( size( M ) ) ;
    
    % RF maps
    for  i = 1 : n.rfmap
      
      % Not selected, empty RF map and continue to next
      if  ~ h.f.UserData.sel( i )
        M{ i } = [ ] ;
        continue
      end
      
      % Locate axes with this map index
      h.ax = findobj( h.f , 'UserData' , i ) ;
      
      % Find baseline data
      h.b = findobj( h.ax , 'Tag' , 'baseline' ) ;
      
      % Initialise logical index matrix to locate baseline values
      B = false( size( M{ i } ) ) ;
      
      % Baseline points
      for  j = 1 : numel( h.b.XData )
        
        % Mark point
        B(  h.b.XData( j ) == coord( i ).x  , ...
            h.b.YData( j ) == coord( i ).y  ) = 1 ;
        
      end % baseline points
      
      % Subtract median baseline
      M{ i } = M{ i }  -  median( M{ i }( B ) ) ;
      
      % Lay down a grid of x- and y-coordinates for each point in M
      X{ i } = repmat( coord( i ).x' , 1 , n.coord.y( i ) ) ;
      Y{ i } = repmat( coord( i ).y  , n.coord.x( i ) , 1 ) ;
      
      % Selected centre falls off the edge of the map
      if  x < min( xlim( h.ax ) )  ||  max( xlim( h.ax ) ) < x  ||  ...
          y < min( ylim( h.ax ) )  ||  max( ylim( h.ax ) ) < y
        
        % Flag as unavailable
        mval( i ) = NaN ;
        
      else
        
        % Find points within 1 standard dev radius of starting centre
        I = sqrt( ( X{ i } - x ) .^ 2 + ( Y{ i } - y ) .^ 2 )  <=  w ;
      
        % Locate the maximum squared value
        mval( i ) = max( M{ i }( I ) .^ 2 ) ;

        % Square root to get magnitude
        mval( i ) = sqrt( mval( i ) ) ;
        
        % Is this zero? Then set to NaN.
        if  mval( i ) == 0  ,  mval( i ) = NaN ;   end
        
      end % is centre on RF map?
      
      % Find valid data points, no NaN or selected baseline values
      I = ~ isnan( M{ i } )  &  ~ B ;
      
      % Keep valid points for least-squares fit
      M{ i } = M{ i }( I )' ;
      X{ i } = X{ i }( I )' ;
      Y{ i } = Y{ i }( I )' ;
      
    end % RF maps
    
    % All mval are NaN or zero
    if  all( isnan( mval ) | mval == 0 )
      
      % Well then set to ones so that the divisions don't produce insane
      % values
      mval( : ) = 1 ;
      
    else
      
      % Any RF map without a normalisation value will get the maximum value
      % found on another map
      mval( isnan( mval ) ) = max( mval ) ;
    
    end % deal with mval
    
    % RF maps, again
    for  i = 1 : n.rfmap
      
      % Not selected, continue to next
      if  ~ h.f.UserData.sel( i )  ,  continue  ,  end
      
      % Normalise RF map
      M{ i } = M{ i }  ./  mval( i ) ;
      
    end % rf maps again
    
    % Collapse data into vectors
    M = [ M{ : } ]' ;
    X = [ X{ : } ]' ;
    Y = [ Y{ : } ]' ;
    
    % Locate maximum squared value
    [ ~ , i ] = max( M .^ 2 ) ;
    
      % Return original for starting amplitude
      a = M( i ) ;
    
    % Define starting coefficients, baseline is zero
    C0 = [ x , y , w , a , 0 ] ;
    
    % Define lower and upper bounds for coefficients
    lb = zeros( 1 , 5 ) ;
    ub = zeros( 1 , 5 ) ;
    
    % Max width will be twice the longest width of the RF map
    ub( 3 ) = 2 * max(  max( X ) - min( X )  ,  ...
                        max( Y ) - min( Y )  ) ;
    
    % Bound Gaussian centre to be within the same distance of zero
    lb( 1 : 2 ) = - ub( 3 ) ;
    ub( 1 : 2 ) = + ub( 3 ) ;
    
    % Minimum and maximum amplitude are bounded from -1 to +1. Same again
    % for baseline.
    lb( 4 : 5 ) = -1 ;  ub( 4 : 5 ) = +1 ;
    
    % Find best fitting coefficients
    C( eind , : ) = lsqcurvefit( @fgauss , C0 , [ X , Y ] , M , ...
      lb , ub , CONST.lop ) ;
    
    % Evaluate best-fitting Gaussian at empirical data points
    G = fgauss( C( eind , : ) , [ X , Y ] ) ;
    
    % Residual sum of squares
    SSres = sum( ( M - G ) .^ 2 ) ;

    % Total sum of squares, remember that starting baseline is the mean
    SStot = sum( ( M - mean( M ) ) .^ 2 ) ;

    % Coefficient of determination
    R2( eind ) = 1  -  SSres ./ SStot ;
    
    % Extract centre and width coefficients
    x = C( eind , 1 ) ;  y = C( eind , 2 ) ;  w = C( eind , 3 ) ;
    
    % Find centre and width markers
    h.c = findobj( h.f , 'Tag' , 'centre' ) ;
    h.w = findobj( h.f , 'Tag' ,  'width' ) ;
    
    % Update to show best fitting choice
    set( h.c , 'XData' , x , 'YData' , y )
    set( h.w , 'Position' , [ x - w , y - w , 2 * [ w , w ] ] )
    
    % Make a new figure to show best-fitting Gaussians. Closing this will
    % unblock program execution.
    h.g = figure( CONST.fig{ : } , 'Tag' , 'gaussproof' , ...
      'DeleteFcn' , @( ~ , ~ ) uiresume( h.f ) , ...
        'KeyPressFcn' , @gkeypress , ...
          'OuterPosition' , h.f.OuterPosition , ...
            'Name' , sprintf( 'R^2: %0.4f' , R2( eind ) ) ) ;
      
    % Shift figure over to the right
    h.g.OuterPosition( 1 ) = sum( h.g.OuterPosition( [ 1 , 3 ] ) ) ;
    
    % Guarantee that figure never goes more than halfway off the right edge
    % of the screen
    h.g.OuterPosition( 1 ) = min(  h.g.OuterPosition( 1 )  ,  ...
      CONST.rgo.ScreenSize( 3 ) - h.g.OuterPosition( 3 ) ./ 2 ) ;
    
    % RF maps
    for  i = 1 : n.rfmap
      
      % Not selected, continue to next
      if  ~ h.f.UserData.sel( i )  ,  continue  ,  end
      
      % Locate axes with this map index
      h.ax = findobj( h.f , 'UserData' , i ) ;
      
      % Make matching axes in the Gaussian figure
      h.gx = subplot( n.ax.row , n.ax.col , i , CONST.ax{ : } , ...
        'YDir' , 'Normal' , 'XLim' , xlim( h.ax ) , ...
          'YLim' , ylim( h.ax ) ) ;
      
      axis equal tight
      colorbar
      
      % Turn off the AxesToolbar object
      h.gx.Toolbar.Visible = 'off' ;
      
      % Define vector of x and y marginal coordinates that span the same
      % range as the original RF map at high resolution
      X = xlim( h.ax ) ;  X = ( X( 1 ) : 0.01 : X( 2 ) )' ;
      Y = ylim( h.ax ) ;  Y =   Y( 1 ) : 0.01 : Y( 2 )    ;
      
      % Number of points
      n.gaus.x = numel( X ) ;
      n.gaus.y = numel( Y ) ;
      
      % Copy for each point on Gaussian curve
      X = repmat( X , 1 , n.gaus.y ) ;
      Y = repmat( Y , n.gaus.x , 1 ) ;
      
      % Evaluate best-fitting gaussian at all points
      G = fgauss(  C( eind , : )  ,  [ X( : ) , Y( : ) ]  ) ;
      
      % Reshape into a matrix
      G = reshape( G , size( X ) ) ;
      
      % Plot best-fitting gaussian
      imagesc( h.gx , X( : , 1 ) , Y( 1 , : ) ,G' , ...
          'PickableParts' , 'none' )
      
    end % rf maps
    
    % Wait for user to approve or reject the fit
    uiwait( h.f )
    
    % Check whether figure is still live
    if  ~ ishandle( h.f )
      
      % Make sure proofing figure is gone
      if  ishandle( h.g )
        h.g.DeleteFcn = [ ] ;
        close( h.g )
      end
      
      % Return empty
      C = [] ;  E = [] ;  R2 = [] ;  coord = [] ;  SEL = [] ;  MAP = [] ;
      return
    
    end % main fig closed
    
    % New electrode index given
    if  h.f.UserData.new_eind
      eind = h.f.UserData.new_eind ; % Copy to master eind variable
      h.f.UserData.new_eind = 0 ;    % Lower flag
      delete( h.f.Children )         % Kill axes
      continue                       % New iteration of electrode loop
    end
    
    % Approved, increment electrode index
    if  h.f.UserData.appflg  ,  eind = eind  +  1 ;  end
    
    % Kill axes
    delete( h.f.Children )
    
  end % electrode loop
  
  % Done, kill main figure
  close( h.f )
  
end % makrfgauss


%%% Sub-routines %%%

% 2D isotropic Gaussian function, compatible with 'fun' input argument of
% lsqcurvefit. Returns values y given coefficient vector c and x is the N
% by 2 array of visual coordinates, in degrees. c contains [ x , y , w , 
% a , b ], the x and y centre coordinate, width as standard dev, amplitude,
% and baseline.
function  y = fgauss( c , x )
  
  % Calculate the difference between data and central coordinate , bsx
  % operation
  z = x  -  c( [ 1 , 2 ] ) ;
  
  % Inverse variance
  iv = 1 ./ c( 3 ) ^ 2 ;
  
  % Inverse covariance matrix
  C = [ iv , 0 ; 0 , iv ] ;
  
  % Evaluate function
  y = c( 4 ) * exp( -0.5 * dot( z * C , z , 2 ) )  +  c( 5 ) ;
  
end % fgauss


% Calculate empirical firing rate of multi-unit spikes from electrode eid
% as function of mapping stimulus coordinate, for each RF mapping session
% in d. Returns 1 x n.rfmap cell array M, containing a matrix of average
% firing rates. M{ i }( xi , yi ) is the mean firing rate in session d( i )
% at coordinate coord( i ).x( xi ) and coord( i ).y( yi ). NaN values
% returned at coordinates with no trials.
function  M = calcrfmap( VARNAM , coord , d , eid )
  
  % Number of RF maps
  n.map = numel( d ) ;
  
  % Shorthand for x- and y-coordinate variable names
  xvar = VARNAM{ 1 } ;
  yvar = VARNAM{ 2 } ;
  
  % Allocate output
  M = cell( 1 , n.map ) ;
  
  % RF maps
  for  i = 1 : n.map
    
    % Go to next mapping session if this electrode not supported
    if  ~ any( d( i ).electrodes == eid )  ,  continue  ,  end
    
    % Number of marginal x and y coordinates
    n.x = numel( coord( i ).x ) ;
    n.y = numel( coord( i ).y ) ;
    
    % Allocate map
    M{ i } = zeros( n.x , n.y ) ;
    
    % Allocate trial count
    n.t = zeros( size( M{ i } ) ) ;
    
    % Trials
    for  t = 1 : d( i ).numtrials
      
      % Find indices for x- and y-coordinate on this trial
      xi = d( i ).var.( xvar )( t )  ==  coord( i ).x ;
      yi = d( i ).var.( yvar )( t )  ==  coord( i ).y ;
      
      % Locate all spikes from this electrode on this trial
      spk = d( i ).spike.electrode{ t }  ==  eid ;
      
      % Accumulate spike count
      M{ i }( xi , yi ) = M{ i }( xi , yi )  +  sum( spk ) ;
      
      % Count trial
      n.t( xi , yi ) = n.t( xi , yi )  +  1 ;
      
    end % trials
    
    % Find coordinates with trials
    j = 0  <  n.t ;
    
    % Compute average firing rate
    M{ i }( j ) = M{ i }( j )  ./  n.t( j ) ;
    
    % Set coords with no trials to NaN
    M{ i }( ~ j ) = NaN ;
    
    % If the entire RF map is zero then no spikes were counted, return an
    % empty matrix because this data is not valid
    if  all(  ~ isfinite( M{ i }( : ) )  |  M{ i }( : ) == 0  )
      M{ i } = [ ] ;
    end
    
  end % RF maps
  
end % calcrfmap


% Update figure name to show where we're at. Returns true if figure was
% closed.
function  ish = fignam( f , i , E , eind , eid , msg )
  
  % Check if handle is still live, quit if not
  ish = ~ ishandle( f ) ;
  if  ish  ,  return  ,  end
  
  % Update figure
  f.UserData.step = i ;
  f.Name = sprintf( 'EID %d (%d of %d) Step %d) %s' , E( eind ) , eind ,...
    eid , f.UserData.step , msg ) ;
  
end % fignam


%%% Callbacks %%%

% KeyPress callback for main figure. 
function  fkeypress( f , k )
  
  % Responds to <Enter> i.e. <Return> or <Escape> keys
  switch  k.Key
    
    % Unblock program execution
    case  'return'
      
      % Respond according to current step
      switch  f.UserData.step
        
        % No special action, simply unbolock
        case  1  ,  uiresume( f )
          
        % Check that baseline values were chosen from all selected RF maps
        case  2  ,  sc = findobj( f , 'Tag' , 'baseline' ) ;
          
          % If no scatter objects are empty then unblock
          if  ~ any( cellfun( @isempty , { sc.XData } ) )
            
            uiresume( f )
            
          % Tell the user what's up
          else
            msgbox( 'Please select baseline values for all RF maps' , ...
              'modal' )
          end
          
        % Make sure that a centre was selected
        case  3  ,  ln = findobj( f , 'Tag' , 'centre' ) ;
          
          % Point selected
          if  ~ isempty( ln( 1 ).XData )
            
            uiresume( f )
            
          % Tell the user what's up
          else
            msgbox( 'Please select starting centre for Gaussian' , ...
              'modal' )
          end
          
        % Make sure that a width was set
        case  4  ,  rt = findobj( f , 'Tag' , 'width' ) ;
          
          % Width selected
          if  ~ all( rt( 1 ).Position == 0 )
            
            uiresume( f )
            
          % Tell the user what's up
          else
            msgbox( 'Please set starting width for Gaussian' , ...
              'modal' )
          end
          
        % User accepts fit, raise flag and close Gaussian proof figure
        case  5  ,  f.UserData.appflg = true ;
                    close( findobj( 'Tag' , 'gaussproof' ) ) ;
                    uiresume( f )
          
      end % respond by step number
      
    % Reset current step
    case  'escape'
      
      % Shift key is down
      if  any( strcmp( k.Modifier , 'shift' ) )
        
        % Get input
        eind = inputdlg( 'Please enter electrode index.' ) ;
        
        % Convert to double
        eind = str2double( eind ) ;
        
        % Invalid input
        if  ~ isscalar( eind )  ||  ~ isfinite( eind )  ||  eind < 0  ||...
            eind > f.UserData.maxein  ||  mod( eind , 1 )
          
          errordlg( 'Not a valid electrode index.' )
          
        % Valid input
        else
          
          % Raise new eind flag
          f.UserData.new_eind = eind ;
          
          % We must unblock execution to go to this electrode
          uiresume( f )
          
        end % validate input
        
        % Done
        return
        
      end % Shift-Escape
      
      % Respond according to current step
      switch  f.UserData.step
        
        % De-select everything
        case  1  ,  f.UserData.sel( : ) = 0 ;
                    arrayfun( @( a ) title( a , '' ) , f.Children ) ;
          
        case  2  ,  sc = findobj( f , 'Tag' , 'baseline' ) ;
                    set( sc , 'XData' , [ ] , 'YData' , [ ] )
                    
        case  3  ,  ln = findobj( f , 'Tag' , 'centre' ) ;
                    set( ln , 'XData' , [ ] , 'YData' , [ ] )
                    
        case  4  ,  rt = findobj( f , 'Tag' , 'width' ) ;
                    set( rt , 'Position' , zeros( 1 , 4 ) )
        
        % Guarantee approval flag is low, kill proof figure, resume
        % execution
        case  5
          
          % Check that the user really means this
          switch  questdlg( 'Would you like to run fit again?' )
            case  'Yes'  ,  f.UserData.appflg = false ;
                            close( findobj( 'Tag' , 'gaussproof' ) ) ;
                            uiresume( f )
          end
          
      end % respond by step number
      
  end % respond
  
  % Store name of latest key press
  f.UserData.key = k.Key ;
  
end % figure keypress


% Key press callback for Gaussian proofing figure
function  gkeypress( f , k )
  
  % Map key pressed to action
  switch  k.Key
    
    % Hitting <Enter> is the same as accepting valitation. Close figure to
    % trigger approval of the fit.
    case  'return'  ,  close( f )
      
    % Escape key, pass this on to the main fig and activate its own key
    % press callback
    case  'escape'
      
      m = findobj( 'Tag' , 'mainfig' ) ;
      m.KeyPressFcn( m , k )
      
  end % map key to action
  
end % gkeypress


% Step 1 mouse click callback for axes objects.
function  a1click( ax , ~ )
  
  % Point to figure
  f = ax.Parent ;
  
  % Valid RF maps
  v = f.UserData.val ;
  
  % RF map index for this axes
  i = ax.UserData ;
  
  % Look for state of shift key. True means down.
  if  strcmp( f.SelectionType , 'extend' )
    
    % Find axes in figure
    ax = findobj( f , 'Type' , 'Axes' ) ;
    
    % If all valid axes are selected then un-select all of them
    if  all( f.UserData.sel( v ) )
      
      f.UserData.sel( v ) = 0 ;
      arrayfun( @( a ) title( a , '' ) , ax ) ;
      
    % Otherwise, select all
    else
      
      f.UserData.sel( v ) = 1 ;
      arrayfun( @( a ) title( a , 'SELECTED' ) , ax ) ;
      
    end
    
  % Shift is up and this axes if valid
  elseif  v( i )
    
    % Reverse selection state
    f.UserData.sel( i ) = ~ f.UserData.sel( i ) ;
    
    % Swap title
    if  f.UserData.sel( i )
      title( ax , 'SELECTED' )
    else
      title( ax , '' )
    end
    
  end % shift is down
  
end % a1click


% Step 2 mouse click callback for axes objects.
function  a2click( ax , coord , M )
  
  % Parent figure
  f = ax.Parent ;
  
  % Gather points according to the type of click
  switch  f.SelectionType
    
    % Right- or control-click to initiate sampling
    case  'alt'
      
      % Set this axes as current
      axes( ax )
      
      % Sample two points
      [ x , y ] = ginput( 2 ) ;
      
      % Two points not given, abort
      if  numel( x ) ~= 2  ,  return  ,  end
      
      % Find points within selected bounds
      xi = min( x ) <= coord.x  &  coord.x <= max( x ) ;
      yi = min( y ) <= coord.y  &  coord.y <= max( y ) ;
      
      % Get rid of NaN values, bsx
      i = xi'  &  yi  &  isfinite( M ) ;
      
      % No valid data
      if  ~ any( i( : ) )  ,  return  ,  end
      
      % Get subscript indices
      [ xi , yi ] = ind2sub( size( M ) , find( i )' ) ;
    
    % Anything else, we take the click as a plain left-click to select a
    % single point
    otherwise
  
      % Get the x- and y-coordinate of the click
      x = ax.CurrentPoint( 1 , 1 ) ;
      y = ax.CurrentPoint( 1 , 2 ) ;

      % Find nearest point on the map
      [ ~ , xi ] = min( abs( coord.x - x ) ) ;
      [ ~ , yi ] = min( abs( coord.y - y ) ) ;
      
      % This point has no data, return with no changes
      if  isnan( M( xi , yi ) )  ,  return  ,  end
  
  end % get coordinates
  
  % Find scatter object
  sc = findobj( ax , 'Tag' , 'baseline' ) ;
  
  % Number of selected points
  np = numel( xi ) ;
  
  % Points to keep, raise flag if kept
  k = false( 1 , np ) ;
  
  % Selected points
  for  i = 1 : np
  
    % Fetch selected point in RF map
    x = coord.x( xi( i ) ) ;
    y = coord.y( yi( i ) ) ;

    % This point has not been added already, raise flag
    k( i ) = ~ any( sc.XData == x  &  sc.YData == y ) ;
  
  end % selected points
  
  % No points kept
  if  ~ any( k )  ,  return  ,  end
  
  % Get new points
  x = coord.x( xi( k ) ) ;
  y = coord.y( yi( k ) ) ;
    
  % Add new points
  set( sc , 'XData' , [ sc.XData , x ] , 'YData' , [ sc.YData , y ] )
  
end % a2click


% Step 3 mouse click callback sets centre coordinate
function  a3click( ax , ~ )
  
  % Find all centre objects
  h = findobj( ax.Parent , 'Tag' , 'centre' ) ;
  
  % Get the x- and y-coordinate of the click
  x = ax.CurrentPoint( 1 , 1 ) ;
  y = ax.CurrentPoint( 1 , 2 ) ;
  
  % Set the centre of all line objects
  set( h , 'XData' , x , 'YData' , y )
  
end % a3click


% Step 4 mouse click callback sets gaussian width parameter
function  a4click( ax , ~ )
  
  % Find centre marker
  c = findobj( ax , 'Tag' , 'centre' ) ;
  
  % Width circles
  w = findobj( ax.Parent , 'Tag' , 'width'  ) ;
  
  % Get the x- and y-coordinate of the click
  x = ax.CurrentPoint( 1 , 1 ) ;
  y = ax.CurrentPoint( 1 , 2 ) ;
  
  % Compute standard deviation
  sd = sqrt(  ( c.XData - x ) ^ 2  +  ( c.YData - y ) ^ 2  ) ;
  
  % Set rectangle position
  set( w , 'Position', [ c.XData - sd , c.YData - sd , 2 * [ sd , sd ] ] );
  
end % a4click

