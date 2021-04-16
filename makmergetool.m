
classdef  makmergetool  <  handle
% 
% obj = makmergetool
% 
% MET Analysis Kit, pre-processing. Returns the object handle to a GUI tool
% that allows a set of spike clusters to be merged together, those from one
% electrode, for example.
% 
% 
% GUI controls:
% 
%   Select a cluster by clicking its average waveform, any of it's points
%   in 3D PCA space, or one of the density histograms in 2D PCA space. Two
%   clusters can be selected at once, one after another. All clusters can
%   be de-selected by clicking the background of any axes. Click the figure
%   background to toggle the rotation tool on and off ; right-click or
%   control-click the figure background to toggle rotation and de-select
%   all clusters ; shift-click the background to toggle 2D density
%   histogram y-axis between waveform count (i.e. number of waveforms per
%   bin) and probability-density (i.e. area under histogram sums to 1).
%   Once two clusters are selected, click the Merge button to join them
%   together into a single cluster. If one cluster is selected, click
%   Reject to discard all of its waveforms. Click the Auto button when no
%   cluster is selected, or if two clusters are selected, to automatically
%   choose the next pair of clusters ; the pair with the highest connection
%   strength is chosen. If one cluster is selected, click the Auto button
%   to automatically find the cluster with the highest connection
%   strength to the first one. The Reset button causes the merge method to
%   request that the calling function re-run merge with cluster data
%   following automated merging in makprep, while the Init button requests
%   inition cluster information prior to any merging. Cutoff causes the
%   merge method to request that automated merging is performed from the
%   intitial cluster set up to the given connection strength cutoff value
%   in the neighbouring edit box. Click the Done button to finish manual
%   merging
% 
% 
% There is a single method:
%
%   [ man , E , c ] = h.merge(  t  ,  w  ,  p  ,  E  ,  c  )
% 
% that initialises the tool using a new set of spikes with spike times t in
% seconds, aligned waveforms w, principal components p, interface energy
% matrix E, and spike cluster assignments c ; t and c are vectors with one
% element per spike, w and p columns are indexed by spike, and E is a
% square matrix. Such data is obtained from the d and s output data structs
% of makprep. Hence, the tool is set using the results of automated spike
% clustering from an initial set of clusters. The function blocks until
% manual clustering is complete, at which point it returns a list of manual
% cluster mergers in man, the final interface energy in E, and the final
% spike cluster assignments in c ; man is a 2 x N matrix with columns
% indexed by merger in the order of occurrence, and rows containing the two
% indices of the merged clusters. During a merger, the high-numbered
% cluster's spikes are added to the low-numbered cluster. If clustering is
% aborted by the user then the output arguments are all empty matrices.
% 
% Written by Jackson Smith - February 2018 - DPAG , University of Oxford
% 
  
  
  % Constant object properties
  properties  (  Constant  )
    
    % Bins used to find ISI histograms , in seconds
    x = 0 : 0.001 : 0.1 ;
    
    % Unselected marker type
    umk = '.' ;
    
    % Unselected marker size
    ums = 18 ;
    
    % Unselected waveform line width
    ulw = 1 ;
    
    % Selected marker type
    smk = 'o' ;
    
    % Selected marker size
    sms =  4 ;
    
    % Selected waveform line width
    slw = 2.5 ;
    
    % PCA line width
    plw = 1.5 ;
    
    % Histogram normalisation modes
    nor = { 'pdf' , 'count' } ;
    
    % Patch face alpha , default
    palpha = 0.5 ;
    
    % Patch face alpha , unselected
    upa = 0.4 ;
    
    % Patch face alpha , selected
    spa = 0.6 ;
    
  end % constant properties


  % Variable object properties
  properties
    
    % Figure handle
    f
    
    % Axes handles for all average waveforms
    aavg
    
    % Principal component axes
    apca
    
    % Selected average waveforms axes
    asel
    
    % ISI histogram axes
    aisi
    
    % All PCA 2D histogram axes
    aph2
    
    % Selected PCA 2D histogram axes
    ash2
    
    % Current 2D density histogram normalisation type
    inor = 1 ;
    
    % Title text
    txttit
    
    % Button for merging pair of cluster
    bmerge
    
    % Button rejects selected cluster
    breject
    
    % Automatically select the next pair of clusters to merge
    bauto
    
    % Button resets clusters to their original spike assignments following
    % automated merging
    breset
    
    % Button resets clusters all the way back to the initial cluster set,
    % for manual merging of all clusters
    binit
    
    % Button for applying connection strength cutoff value provided in the
    % accompanying edit box
    bcutoff
    ecutoff
    
    % Button finishes manual cluster merging
    bdone
    
    % Title header
    thead
    
    % Average waveform handles per cluster , with standard deviation
    % patches
    avg
    
    % PCA clouds per cluster
    pca
    
    % Average waveforms of selected clusters , with standard deviation
    % patches
    sel
    
    % Inter-spike-interval histograms of selected clusters
    isi
    
    % PCA 2D histograms
    ph2
    
    % Selected PCA 2D histograms
    sh2
    
    % The number of clusters
    n
    
    % Unique cluster indices
    u
    
    % Spike times
    t
    
    % Waveforms
    w
    
    % Principal components
    p
    
    % Interface energy matrix
    E
    
    % Spikes per cluster
    spc
    
    % Spike cluster id's
    c
    
    % Connection strength
    J
    
    % Low-index selected cluster
    c1
    
    % High-index selected cluster
    c2
    
    % Merge list
    man
    
    % Return flag , this is a string initialised to empty. If not empty
    % then it is returned as the first argument from merge
    retflg
    
  end % variable properties
  
  
  % Object methods
  methods
    
    function  h = makmergetool
    %
    % Create object
    %
      
      % Get a basic figure without the menu bar or dock controls , keep
      % invisible until set up
      h.f = figure (  'DockControls' ,  'off' ,  'MenuBar' ,  'none' ,  ...
        'ToolBar' ,  'figure' ,  'Visible' ,  'off'  ,  ...
          'Name'  ,  'MAK spike cluster merge tool'  ,  ...
            'NumberTitle'  ,  'off'  ,  'ButtonDownFcn'  ,  @figure_cb  ) ;
      
      % Fill the screen
      h.f.OuterPosition = get (  groot  ,  'ScreenSize'  ) ;
      
      % Set rotation object callback
      robj = rotate3d (  h.f  ) ;
      robj.ButtonDownFilter = @robj_cb ;
      
      % Create global average waveform axes
      h.aavg = subplot (  2  ,  3  ,  1  ,  'Parent'  ,  h.f  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
        
      % 3D PCA scatter plot axes
      h.apca = subplot (  2  ,  3  ,  2  ,  'Parent'  ,  h.f  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
        
      % Global 2D PCA density histogram plot axes
      h.aph2 = subplot (  2  ,  3  ,  3  ,  'Parent'  ,  h.f  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
      
      % Selected average waveform axes
      h.asel = subplot (  2  ,  3  ,  4  ,  'Parent'  ,  h.f  ,  ...
        'ActivePositionProperty'  ,  'position'  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
      
      % ISI histogram axes
      h.aisi = subplot (  2  ,  3  ,  5  ,  'Parent'  ,  h.f  ,  ...
        'ActivePositionProperty'  ,  'position'  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
      
      % Selected 2D PCA density histogram plot axes
      h.ash2 = subplot (  2  ,  3  ,  6  ,  'Parent'  ,  h.f  ,  ...
        'ActivePositionProperty'  ,  'position'  ,  ...
          'ButtonDownFcn'  ,  @axes_cb  ) ;
      
      % Set certain common properties
      a = [ h.aavg , h.apca , h.asel , h.aisi , h.aph2 , h.ash2 ] ;
      set (  a  ,  'TickDir'  ,  'out'  ,  'LineWidth'  ,  1  ,  ...
        'FontSize'  ,  12  )
      for  a = a
        grid (  a  ,  'on'  )
        hold (  a  ,  'on'  )
         box (  a  ,  'on'  )
      end
      
      % Axis labels
      xlabel (  h.aavg  ,  'Samples'  )
      xlabel (  h.apca  ,  'PCA 1'  )
      xlabel (  h.aph2  ,  'PCA 1'  )
      xlabel (  h.asel  ,  'Samples'  )
      xlabel (  h.aisi  ,  'Inter-spike-interval ( sec )'  )
      xlabel (  h.ash2  ,  'PCA 1'  )
      
      ylabel (  h.aavg  ,  'Average waveform (micro-volt)'  )
      ylabel (  h.apca  ,  'PCA 2'  )
      ylabel (  h.aph2  ,  'PCA 2'  )
      ylabel (  h.asel  ,  'Avg wave ( uv )'  )
      ylabel (  h.aisi  ,  h.nor{ h.inor }  )
      ylabel (  h.ash2  ,  'PCA 2'  )
      
      zlabel (  h.apca  ,  'PCA 3'  )
      zlabel (  h.aph2  ,  h.nor{ h.inor }  )
      zlabel (  h.ash2  ,  h.nor{ h.inor }  )
      
      % Special properties
      axis (  h.apca  ,  'equal'  )
      
      % Create buttons
      h.bmerge = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Enable'  ,  'off'  ,  'Callback'  ,  @bmerge_cb  ,  ...
        'TooltipString'  ,  'Merge selected spike clusters'  ,  ...
        'String'  ,  'Merge'  ,  'Units'  ,  'normalized'  ) ;
      h.breject = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Enable'  ,  'off'  ,  'Callback'  ,  @breject_cb  ,  ...
        'TooltipString'  ,  'Reject selected spike cluster'  ,  ...
        'String'  ,  'Reject'  ,  'Units'  ,  'normalized'  ) ;
      h.bauto = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Callback'  ,  @bauto_cb  ,  ...
        'TooltipString'  ,  sprintf ( ...
          [ 'Automatically select next cluster pair.\n' ...
            'Searches all pairs if zero or two clusters selected.\n' , ...
            'If one cluster selected then it finds the best match.' ] ),...
        'String'  ,  'Auto'  ,  'Units'  ,  'normalized'  ) ;
      h.breset = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Callback'  ,  @breset_cb  ,  'UserData'  ,  'reset'  ,  ...
        'TooltipString'  ,  'Reset to automated clusters assignments'  ,...
        'String'  ,  'Reset'  ,  'Units'  ,  'normalized'  ) ;
      h.binit = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Callback'  ,  @breset_cb  ,  'UserData'  ,  'init'  ,  ...
        'TooltipString'  ,  'Reset to initial clusters assignments'  ,...
        'String'  ,  'Init'  ,  'Units'  ,  'normalized'  ) ;
      h.bcutoff = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Callback'  ,  @breset_cb  ,  'UserData'  ,  'cutoff'  ,  ...
        'TooltipString' , 'Re-merge clusters up to connection strength',...
        'String'  ,  'Cutoff'  ,  'Units'  ,  'normalized'  ) ;
      h.ecutoff = uicontrol (  h.f  ,  'Style'  ,  'edit'  ,  ...
        'Callback'  ,  @ecutoff_cb  ,  'Units'  ,  'normalized'  ,  ...
          'TooltipString'  ,  'Enter new connection strength'  ,  ...
            'String'  ,  '0.1000'  ) ;
      h.bdone = uicontrol (  h.f  ,  'Style'  ,  'pushbutton'  ,  ...
        'Enable'  ,  'off'  ,  'Callback'  ,  @bdone_cb  ,  ...
        'TooltipString'  ,  'Finished merging clusters'  ,  ...
        'String'  ,  'Done'  ,  'Units'  ,  'normalized'  ) ;
      
      % Cutoff edit requires a struct in UserData
      h.ecutoff.UserData = struct (  'str'  ,  h.ecutoff.String  ,  ...
        'num'  ,  str2double( h.ecutoff.String )  ) ;
      
      % Resize buttons according to width of string
      for  a = [ h.bmerge , h.breject , h.bauto , h.breset , h.binit , ...
          h.bcutoff , h.ecutoff , h.bdone ]
        
        a.Position( 3 ) = 1.1  *  a.Extent( 3 ) ;
        
      end % resize buttons
      
      % Align Merge button at top-left of global average waveform axes
      p = h.aavg.Position ;
      h.bmerge.Position( 1 ) = p( 1 ) ;
      h.bmerge.Position( 2 ) = sum ( p(  [ 2 , 4 ]  ) )  -  ...
        h.bmerge.Position( 4 ) ;
      
      % Raise remaining buttons to same vertical position
      p = h.bmerge.Position ;
      for  a = [ h.breject , h.bauto , h.breset , h.binit , h.bcutoff , ...
          h.ecutoff , h.bdone ]
        a.Position( 2 ) = p( 2 ) ;
      end
      
      % Place buttons in sequence from left to right
      for  a = {  [ h.bmerge , h.breject ] ,  [ h.breject , h.bauto ]  ,...
          [ h.bauto , h.breset ] ,  [ h.breset , h.binit ] ,  ...
          [ h.binit , h.bcutoff ] , [ h.bcutoff , h.ecutoff ] , ...
          [ h.ecutoff , h.bdone ]  }
        
        p = a{ 1 }( 1 ).Position ;
        a{ 1 }( 2 ).Position( 1 ) = sum ( p(  [ 1 , 3 ]  ) ) ; %#ok
        
      end % left to right
      
      % Cut height of average waveform axes
      h.aavg.Position( 4 ) = 0.75  *  h.aavg.Position( 4 ) ;
      
      % Make title text object
      p = h.bmerge.Position ;
      h.txttit = uicontrol (  h.f  ,  'Style'  ,  'text'  ,  ...
        'String'  ,  ''  ,  'Units'  ,  'normalized'  ,  ...
          'FontSize' ,  11 ,  'FontWeight' ,  'bold' ,  ...
            'HorizontalAlignment' ,  'left'  ) ;
      h.txttit.Position( 1 ) = p( 1 ) ;
      h.txttit.Position( 2 ) = sum( p( [ 2 , 4 ] ) ) ;
      
      
      % At last , link makmergetool object to its figure
      h.f.UserData = h ;
      
      % Finished , reveal the beauty that is this figure
      h.f.Visible = 'on' ;
      
    end % create object
    
    
    function  delete (  h  )
    %
    % Delete the object , including its figure
    %
      
      delete (  h.f  )
      
    end % delete object
    
    
    function  [ man , Eout , cout , uout ] = merge (  h  ,  ...
        t  ,  w  ,  p  ,  E  ,  spc  ,  c , thead  )
    %
    % Merge function sets up GUI tool using given spike data then blocks
    % until either the Done button is hit or until the figure is closed. If
    % Done is hit then list of manual merges, final energy matrix, spike
    % cluster assignments, and unique cluster numbers are returned. If
    % Reset, Init, or Cutoff are hit then man returns strings 'return',
    % 'init', and 'cutoff' respectively ; for Cutoff, Eout will contain the
    % new connection strength cutoff value.
    %
      
      % Default output
      man = [] ;  Eout = [] ;  cout = [] ;  uout = [] ;
      
      % Purge merger list
      h.man = zeros (  2  ,  0  ) ;
      
      % Initialise return flag
      h.retflg = '' ;
      
      % Store data
      h.t = t ;  h.w = w ;  h.p = p ;  h.E = E ;  h.spc = spc ;  h.c = c ;
      if  7  <  nargin
        h.txttit.String = thead ;
        h.txttit.Position( 3 : 4 ) = 1.1  *  h.txttit.Extent( 3 : 4 ) ;
      else
        h.txttit.String = '' ;
      end
      
      % Check input
      
        % If PCA returned fewer than three components then zero pad missing
        % dimensions
        if  size (  h.p  ,  1  )  <  3
          
          h.p = [  h.p  ;
                   zeros(  3 - size( h.p , 1 )  ,  size( h.p , 2 )  )  ] ;
          
        end % check pca size
      
      % Calculate connection strength
      h.Jrefresh
      
      % Unique cluster indices
      h.u = unique (  c  ) ;
      
      % Number of clusters
      h.n = numel (  h.u  ) ;
      
        % There are no clusters , return empties
        if  ~ h.n  ,  return  ,  end
        
      % Just in case this was missed , delete any trailing graphics objects
      delete (  [ h.avg ; h.pca ; h.sel ; h.isi ; h.ph2 ; h.sh2 ]  )
        
      % Make new graphics object vectors
      h.avg = gobjects (  2 * h.n  ,  1  ) ;
      h.pca = gobjects (  h.n  ,  1  ) ;
      h.sel = gobjects (  2 * h.n  ,  1  ) ;
      h.isi = gobjects (  h.n  ,  1  ) ;
      h.ph2 = gobjects (  h.n  ,  1  ) ;
      h.sh2 = gobjects (  h.n  ,  1  ) ;
      
      % Reset axes colour index
      set (  [ h.aavg , h.apca , h.asel , h.aisi , h.aph2 , h.ash2 ]  , ...
        'ColorOrderIndex'  ,  1  )
        
      % Add cluster graphics objects
      for  i = 1 : h.n
        
        % Find all spikes in this cluster
        j = c  ==  h.u( i ) ;
        
        % Compute average waveform
        wavg = mean (  w( : , j )  ,  2  ) ;
        
        % Compute waveform standard deviation
        wstd =  std (  w( : , j )  ,  0  ,  2  ) ;
        
        % Add line to axes with all average waveforms. Line responds to
        % mouse clicks and knows which cluster it represents.
        h.avg( i ) = plot (  h.aavg ,  wavg ,  'LineWidth' ,  h.ulw ,  ...
          'ButtonDownFcn' ,  @cline_cb ,  'UserData' ,  h.u( i )  ) ;
        
        % Add corresponding standard deviation patch. Make it respond to
        % mouse clicks, knowing which cluster it represents
        h.avg( h.n + i ) = h.getpatch (  h.aavg  ,  wavg  ,  wstd  ,  ...
          h.avg( i ).Color  ,  h.u( i )  ,  h.upa  ) ;
        h.avg( h.n + i ).ButtonDownFcn = @cline_cb ;
        
        % Add average waveform to selected waveform axes. Since nothing is
        % selected to begin with, this is not visible. Knows which cluster
        % it is and selects it following mouse click.
        h.sel( i ) = plot (  h.asel ,  wavg ,  'LineWidth' ,  h.ulw ,  ...
          'UserData' ,  h.u( i ) ,  'Visible' ,  'off'  ,  ...
            'ButtonDownFcn' ,  @cline_cb  ) ;
        
        % Add corresponding standard deviation patch. Make invisible
        h.sel( h.n + i ) = h.getpatch (  h.asel  ,  wavg  ,  wstd  ,  ...
          h.sel( i ).Color  ,  h.u( i )  ,  h.palpha  ) ;
        h.sel( h.n + i ).ButtonDownFcn = @cline_cb ;
        h.sel( h.n + i ).Visible = 'off' ;
        
        % Add line to PCA axes. Again, line responds to mouse clicks and
        % knows which cluster it represents.
        h.pca( i ) = plot3 (  h.apca  ,  ...
          h.p( 1 , j )  ,  h.p( 2 , j )  ,  h.p( 3 , j )  ,  h.umk  ,  ...
            'MarkerSize'  ,  h.ums  ,  'LineWidth'  ,  h.plw  ,  ...
              'ButtonDownFcn'  ,  @cline_cb  ,  'UserData'  ,  h.u( i )  );
            
        % Compute ISI histogram. Since nothing is selected to begin with,
        % this is not visible. Always use default normalisation. Knows its
        % cluster and selects it after mouse click.
        h.isi( i ) = histogram (  h.aisi  ,  diff( h.t( j ) )  ,  h.x  ,...
          'Normalization'  ,  h.nor{ 1 }  ,  'UserData'  ,  h.u( i )  , ...
            'Visible'  ,  'off'  ,  'ButtonDownFcn'  ,  @cline_cb  ) ;
          
        % Compute 2D PCA histogram. Histograms respond to mouse clicks and
        % know which clusters they represent.
        h.ph2( i ) = histogram2 (  h.aph2  ,  h.p( 1 , j )  ,  ...
          h.p( 2 , j )  ,  'Normalization'  ,  h.nor{ h.inor }  ,  ...
            'LineWidth'  ,  h.ulw  ,  'ButtonDownFcn'  ,  @cline_cb  ,  ...
              'UserData'  ,  h.u( i )  )  ;
            
        % Selected 2D PCA histogram. Histograms respond to mouse clicks and
        % know which clusters they represent. Not visible to start with.
        h.sh2( i ) = histogram2 (  h.ash2  ,  h.p( 1 , j )  ,  ...
          h.p( 2 , j )  ,  'Normalization'  ,  h.nor{ h.inor }  ,  ...
            'LineWidth'  ,  h.ulw  ,  'ButtonDownFcn'  ,  @cline_cb  ,  ...
              'UserData'  ,  h.u( i )  ,  'Visible'  ,  'off'  )  ;
        
      end % average waveforms
      
      % Average waveform x-axis limits
      set (  [ h.aavg , h.asel ]  ,  ...
        'xlim'  ,  [ 0 , numel( wavg ) + 1 ]  )
      
      % Re-order st.dev patches and average lines in waveform axes
      h.aavg.Children = [  findobj( h.aavg , 'Type' ,  'line' )  ;
                           findobj( h.aavg , 'Type' , 'patch' )  ] ;
      h.asel.Children = [  findobj( h.asel , 'Type' ,  'line' )  ;
                           findobj( h.asel , 'Type' , 'patch' )  ] ;
      
      % Nothing selected
      h.select( 0 )
      
      % Activate Done button
      h.bdone.Enable = 'on' ;
      
      % Block on GUI tool
      uiwait (  h.f  )
      
      % Figure was not deleted
      if  isgraphics (  h.f  )
      
        % Return flag is empty string
        if  isempty (  h.retflg  )
          
          % Return data
          man = h.man ;
          Eout = h.E ;
          cout = h.c ;
          uout = h.u ;
          
        % Return flag is non-empty string
        else
          
          % So return as first output argument
          man = h.retflg ;
          
          % Returning new connection strength cutoff value
          if  strcmp (  man  ,  'cutoff'  )
            Eout = h.ecutoff.UserData.num ;
          end
          
        end % return flag
      
      end % valid figure handle
      
    end % merge
    
    
    function  [ X , Y ] = getpatchcoords (  ~  ,  wavg  , wstd  )
    %
    % Calculate the x- and y-coordinates for a standard deviation patch
    %
      
      % Number of waveform samples
      N = numel (  wavg  ) ;
      
      % Make x-axis vector
      X = [  1 : N  ,  N : -1 : 1  ]' ;
      
      % Make y-axis vector
      Y = [  wavg + wstd  ; ...
             wavg( end : -1 : 1 )  -  wstd( end : - 1 : 1 )  ] ;
      
    end % getpatchcoords
    
    
    function  p = getpatch (  h  ,  a  ,  wa  ,  ws  ,  C  ,  U  ,  FA  )
    %
    % Returns a patch object that shows the standard deviation around the
    % average waveform.
    %
      
      % Patch coordinates
      [ X , Y ] = h.getpatchcoords(  wa  ,  ws  ) ;
      
      % Build patch object
      p = patch (  X ,  Y ,  C ,  'Parent' ,  a ,  ...
        'FaceAlpha' ,  FA ,  'LineStyle' ,  'none' ,  'UserData' ,  U  ) ;
      
    end % getpatch
    
    
    function  seltit (  h  )
    %
    % Selection title that appears over the global average waveform axes.
    % Expects h.c1  <  h.c2
    %
      
      % Nothing selected
      if  ~ (  h.c1  ||  h.c2  )
        
        tit = sprintf (  'No selection\n'  ) ;
        
      % Both selected
      elseif  h.c1  &&  h.c2
        
        tit = sprintf (  [ 'Cluster 1: %d (%d)\nCluster 2: %d (%d)\n' , ...
          'Connection strength: %f' ]  ,  h.c1  ,  h.spc( h.c1 )  ,  ...
            h.c2  ,  h.spc( h.c2 )  ,  full( h.J( h.c1 , h.c2 ) )  ) ;
          
      % One selected
      elseif  h.c1  &&  ~ h.c2
        
        tit = sprintf (  'Cluster 1: %d (%d)\n\n'  ,  ...
          h.c1  ,  h.spc( h.c1 )  ) ;
        
      % If we get here then there was an error
      else
        
        error (  'MAK:makmergetool:seltit'  ,  ...
          'makmergetool:seltit coding error'  )
        
      end
      
      % Set title
      title (  h.aavg  ,  tit  )
      
    end % seltit
    
    
    function  select (  h  ,  c  )
    %
    % Update cluster selection, adding cluster index c. If c is 0 then
    % selections are cleared.
    %
      
      % Remove highlighting from these clusters
      hrem = [] ;
      
      % Add highlighting to this cluster
      hadd = [] ;
      
      % Update selected cluster indices and get list of which clusters to
      % highlight or un-highlight. First check if we're deselect
      % everything.
      if  c  ==  0
        
        % Get both possible selections
        hrem = [ h.c1 , h.c2 ] ;
        
        % Keep only non-zeros
        hrem( hrem  ==  0 ) = [] ;
        
        % De-select
        h.c1 = 0 ;
        h.c2 = 0 ;
        
        % Disable buttons
        h.bmerge.Enable = 'off' ;
        h.breject.Enable = 'off' ;
        
      % Nothing is selected , select first cluster
      elseif  h.c1  ==  0
        
        % Remember selection & add highlighting
        h.c1 = c ;
        hadd = c ;
        
        % Disable Merge enable reject buttons
        h.bmerge.Enable = 'off' ;
        h.breject.Enable = 'on' ;
        
      % One cluster is selected , select second cluster
      elseif  h.c2  ==  0
        
        % Same cluster selected twice , quit now
        if  h.c1  ==  c  ,  return  ,  end
        
        % Remember selection & add highlighting
        h.c2 = c ;
        hadd = c ;
        
        % Guarantee that c1 is the low-numbered cluster
        c = sort (  [ h.c1 , h.c2 ]  ) ;
        h.c1 = c( 1 ) ;
        h.c2 = c( 2 ) ;
        
        % Disable reject enable merge buttons
        h.bmerge.Enable = 'on' ;
        h.breject.Enable = 'off' ;
        
      % Two clusters are selected , deselect both and select new cluster
      else
        
        % Remove highlighting from these
        hrem = [ h.c1 , h.c2 ] ;
        
        % Add highlighting
        hadd = c ;
        
        % Set cluster selections
        h.c1 = c ;
        h.c2 = 0 ;
        
        % Disable merge and enable reject buttons
        h.bmerge.Enable = 'off' ;
        h.breject.Enable = 'on' ;
        
      end % udate indices , get highlighting sets
      
      % Change highlighting
      h.highlight(  'off'  ,  hrem  )
      h.highlight(   'on'  ,  hadd  )
      
      % Update selection title
      h.seltit
      
    end % select
    
    
    function  highlight (  h  ,  vis  ,  C  )
    %
    % Change highlighting state of selected clusters
    %
      
      % Choose set of line properties
      switch  vis
        case   'on'  ,  mk = h.smk ; ms = h.sms ; lw = h.slw ; pa = h.spa ;
        case  'off'  ,  mk = h.umk ; ms = h.ums ; lw = h.ulw ; pa = h.upa ;
      end
      
      % Change highlighting for each cluster
      for  cind = C
      
        % Find line from global average axes and 2D PCA histograms. Do not
        % include st.dev patches
        l = findobj (  [ h.avg ; h.ph2 ]  ,  'UserData'  ,  cind  ,  ...
          '-not'  ,  'Type'  ,  'patch'  ) ;
        set (  l  ,  'LineWidth'  ,  lw  )
        
        % Now find patch from global average axes and set face alpha
        l = findobj (  h.avg ,  'UserData' ,  cind ,  'Type' ,  'patch'  );
        l.FaceAlpha = pa ;
        
        % Find line from PCA axes
        l = findobj (  h.pca  ,  'UserData'  ,  cind  ) ;
        set (  l  ,  'Marker'  ,  mk  ,  'MarkerSize'  ,  ms  )
        
        % Set visibility of objects from selected waveform, ISI histogram
        % axes, and selected 2D PCA histograms
        l = findobj(  [ h.sel ; h.isi ; h.sh2 ]  ,  'UserData'  ,  cind  );
        set (  l  ,  'Visible'  ,  vis  )
        
      end % clusters
      
    end % highlight
    
    
    function  delclustg (  h  ,  c  )
    %
    % Deletes graphics objects representing spike cluster with index c
    %
      
      % Loop graphics sets
      for  G = { 'avg' , 'pca' , 'sel' , 'isi' , 'ph2' , 'sh2' }
        
        % Get graphics vector name
        g = G{ 1 } ;
        
        % Find the objects representing spike cluster c
        gob = findobj (  h.( g )  ,  'UserData'  ,  c  ) ;
        
        % Get their location in graphics handle vector
        i = false ( size(  h.( g )  ) ) ;
        
        for  j = 1 : numel( gob )
          i = i  |  ( h.( g )  ==  gob( j ) ) ;
        end
        
        % Remove them from the graphics set
        h.( g )( i ) = [] ;
        
        % Delete objects
        delete (  gob  )
        
      end % graphics sets
      
    end % delclustg
    
    
    function  delclustdat (  h  ,  c  )
    %
    % Delete cluster data for cluster with index c
    %

      % Remove cluster from list of cluster indices
      h.u( h.u == c ) = [] ;

      % Decrement number of clusters
      h.n = h.n  -  1 ;

      % Zero number of spikes
      h.spc( c ) = 0 ;

      % Remove cluster from interface energy matrix
      h.E( c ,  : ) = 0 ;
      h.E(  : , c ) = 0 ;
      h.E( c , c ) = 1 ;

      % Recompute connection strength matrix
      h.Jrefresh
      
    end % delclustdat
    
    
    function  update2dhist (  h  ,  a  ,  c  ,  i  )
    % 
    % Updates 2D PCA histogram so that it contains new cluster data with
    % appropriate bin limits
    % 
      
      % Locate histogram , a child of axes a
      l = findobj (  a  ,  'UserData'  ,  h.c1  ) ;
      
      % Set data to new set of spikes for specific PCA components
      l.Data = h.p( c , i )' ;
      
      % Determine new set of bin edges
      [ ~ , l.XBinEdges , l.YBinEdges ] = ...
        histcounts2 (  l.Data( : , 1 )  ,  l.Data( : , 2 )  ) ;
      
    end % update2dhist
    
    
    function  Jrefresh (  h  )
    %
    % Compute connection strength matrix from current interface energy
    % matrix
    %
      
      % Compute new connection strengths
      h.J = makconnstrength (  h.E  ,  h.spc  ) ;
      
      % Find all elements on diagonal and in lower-triangular portion of
      % the matrix , or those with NaN values
      i = tril(  true( size( h.J , 2 ) )  )  |  isnan ( h.J ) ;
      
      % Mask out these elements
      h.J( i ) = 0 ;
      
    end % J_update
    
    
  end % methods
  
  
end % makmergetool


%%% Callbacks %%%

% Figure ButtonDownFcn toggles 3D rotation tool on and off. If user
% right-clicks then all clusters are de-selected as well.
function  figure_cb (  f  ,  ~  )
  
  % Get makmergetool object
  h = f.UserData ;
  
  % Choose behaviour based on mouse-click type i.e. selection type
  switch  f.SelectionType
    
    % Simple left-click , toggle rotation only
    case  'normal'  ,  rotate3d (  f  )
      
    % Right- or control-click , toggle rotation and de-select clusters
    case     'alt'  ,  rotate3d (  f  )  ,  h.select( 0 )
      
    % Shift-click , toggle 2D PCA density histogram normalisation
    case  'extend'
      
      % Toggle normalisation index
      h.inor = mod (  h.inor  ,  numel(  h.nor  )  )  +  1 ;
      
      % Set 2D PCA density histogram normalisation
      set (  [ h.ph2 , h.sh2 ]  ,  'Normalization'  ,  h.nor{ h.inor }  )
      
      % And set z-axis label on corresponding axes
      zlabel (  h.aph2  ,  h.nor{ h.inor }  )
      zlabel (  h.ash2  ,  h.nor{ h.inor }  )
    
  end % selection type
  
end % figure_cb


% Rotation object ButtonDownFilter callback. If figure background was
% clicked then allow its ButtonDwnFcn to execute. In this case, toggle the
% rotation tool on and off.
function  res = robj_cb (  obj  ,  ~  )
  
  % Allow ButtonDownFcn to run if the object is a figure
  res = isa (  obj  ,  'matlab.ui.Figure'  ) ;
  
end % robj_cb


% If user clicks background of the global average waveform or pca axes then
% selection is cleared
function  axes_cb (  a  ,  ~  )

  % Get makmergetool object
  h = a.Parent.UserData ;
  
  % Clear selection
  h.select( 0 )

end % axes_cb


% Merge button
function  bmerge_cb (  b  ,  ~  )
  
  % makmergetool object
  h = b.Parent.UserData ;
  
  % No selection or only one selected , quit now. We should never get here.
  if  h.c1  ==  0  ||  h.c2  ==  0
    warning (  'MERGE button coding error'  )
    return
  end
  
  % Add to merge list
  h.man = [  h.man  ,  [ h.c1 ; h.c2 ]  ] ;
  
  % Find all spikes belonging to high-index cluster
  i = h.c  ==  h.c2 ;
  
  % Assign to low-index cluster
  h.c( i ) = h.c1 ;
  
  % Compute new intra-cluster energy
  h.E( h.c1 , h.c1 ) = ...
    h.E( h.c1 , h.c1 )  +  h.E( h.c2 , h.c2 )  +  h.E( h.c1 , h.c2 ) ;

  % Get linear indices of for inter-cluster energies between the low
  % cluster and all others -- with the exception of the high cluster.
  % Takes an 'L' shape out of the upper triangle of the raw energy table.
  [ i1 , i2 ] = makcind (  h.c1  ,  h.c2  ,  size( h.E , 2 )  ) ;

  % Compute new inter-cluster energies
  h.E( i1 ) = h.E( i1 )  +  h.E( i2 ) ;
  
  % Add spikes to low-index cluster count
  h.spc( h.c1 ) = h.spc( h.c1 )  +  h.spc( h.c2 ) ;
  
  % Delete high-index cluster data , this updates connection strength
  % matrix
  h.delclustdat(  h.c2  )
  
  % Find all spikes on low-index cluster
  i = h.c  ==  h.c1 ;
  
  % Compute average waveform , standard deviation and patch coordinates
  wavg = mean (  h.w( : , i )  ,  2  ) ;
  wstd =  std (  h.w( : , i )  ,  0  ,  2  ) ;
  [ X , Y ] = h.getpatchcoords (  wavg  ,  wstd  ) ;
  
  % Update global axes average waveform
  l = findobj (  h.avg  ,  'UserData'  ,  h.c1  ,  'Type'  ,   'line'  ) ;
  l.YData = wavg ;
  l = findobj (  h.avg  ,  'UserData'  ,  h.c1  ,  'Type'  ,  'patch'  ) ;
  set (  l  ,  'XData'  ,  X  ,  'YData'  ,  Y  )
  
  % Update selected axes average waveform
  l = findobj (  h.sel  ,  'UserData'  ,  h.c1  ,  'Type'  ,   'line'  ) ;
  l.YData = wavg ;
  l = findobj (  h.sel  ,  'UserData'  ,  h.c1  ,  'Type'  ,  'patch'  ) ;
  set (  l  ,  'XData'  ,  X  ,  'YData'  ,  Y  )
  
  % Update PCA line
  l = findobj (  h.pca  ,  'UserData'  ,  h.c1  ) ;
  set (  l  ,  'XData'  ,  h.p( 1 , i )  ,  'YData'  ,  h.p( 2 , i )  , ...
    'ZData'  ,  h.p( 3 , i )  )
  
  % And update ISI histogram
  l = findobj (  h.isi  ,  'UserData'  ,  h.c1  ) ;
  l.Data = diff (  h.t( i )  ) ;
  
  % Update 2D PCA histograms
  h.update2dhist(  h.aph2  ,  [ 1 , 2 ]  ,  i  )
  h.update2dhist(  h.ash2  ,  [ 1 , 2 ]  ,  i  )
  
  % Remember low- and high-index clusters
  c = [ h.c1 , h.c2 ] ;
  
  % Deselect everything
  h.select(  0  )
  
  % Delete graphics objects for high-index cluster
  h.delclustg(  c( 2 )  )
  
  % Reselect low-index cluster
  h.select(  c( 1 )  )
  
end % bmerge_cb


% Reject button
function  breject_cb (  b  ,  ~  )
  
  % makmergetool object
  h = b.Parent.UserData ;
  
  % No selection or two are selected , quit now. We should never get here.
  if  h.c1  ==  0  ||  h.c2  ~=  0
    warning (  'Reject button coding error'  )
    return
  end
  
  % Find all spikes belonging to selected cluster
  i = h.c  ==  h.c1 ;

  % Assign null cluster index
  h.c( i ) = 0 ;
  
  % Remember cluster number
  c = h.c1 ;
  
  % Delete cluster data
  h.delclustdat(  h.c1  )
  
  % Nothing selected
  h.select(  0  )
  
  % Delete graphics objects
  h.delclustg(  c  )
  
end % breject_cb


% Automatically select the next cluster pair with the strongest connection
% strength
function  bauto_cb (  b  ,  ~  )
  
  % makmergetool object
  h = b.Parent.UserData ;
  
  % Quit if there are fewer than two clusters
  if  h.n  <  2  ,  return  ,  end
  
  % How many clusters are currently selected?
  nsel = sum (  0  <  [ h.c1 , h.c2 ]  ) ;
  
  % Find cluster pair. If zero or two clusters are currently selected then
  % search all pairs.
  if  any (  nsel  ==  [ 0 , 2 ]  )
  
    % Get linear index of cluster with highest connection strength
    [ ~ , i ] = max (  h.J( : )  ) ;

    % Get row and column i.e. cluster indices
    [ r , c ] = ind2sub (  size( h.J )  ,  i  ) ;
    
  % If one cluster is selected then find its nearest neighbour
  else
    
    % Gather two connection strengths
    j = [ 0 , 0 ] ;
    
    % Search connection strength matrix along selected cluster's row for a
    % column index and along its column for a row index
    [ j( 1 ) , c ] = max (  h.J( h.c1 , : )  ) ;
    [ j( 2 ) , r ] = max (  h.J( : , h.c1 )  ) ;
    
    % Maximum connection strength was found along selected cluster's row
    if  j( 1 )  >  j( 2 )
      
      % So we keep the column index , and say which row it came from
      r = h.c1 ;
      
    % Max conn. strength found along selected's column
    else
      
      % So we keep the row index , and say which column it came from
      c = h.c1 ;
      
    end % In this way , we keep the index of the cluster with the strongest
        % connection strength to the selected cluster
    
  end % find pair
  
  % Deselect
  h.select(  0  )
  
  % Automatic selection
  h.select(  r  )
  h.select(  c  )
  
end % bauto_cb


% Reset spike cluster assignments by raising return flag with non-empty
% string taken from button's UserData , then execute the Done button
function  breset_cb (  b  ,  ~  )
  
  % makmergetool object
  h = b.Parent.UserData ;
  
  % Raise return flag
  h.retflg = b.UserData ;
  
  % Find Done button
  b = findobj (  h.f  ,  'Type'  ,  'uicontrol'  ,  ...
    'Style'  ,  'pushbutton'  ,  'String'  ,  'Done'  ) ;
  
  % Execute Done button
  b.Callback( b , [] )
  
end % breset_cb


% Check new input to make sure that it is a scalar, positive, real number
% between 0 and 1
function  ecutoff_cb (  e  ,  ~  )
  
  % Convert string to double numeric value
  d = str2double (  e.String  ) ;
  
  % Check value
  if  ~ isscalar (  d  )  ||  isnan (  d  )  ||  isinf (  d  )  ||  ...
      ~ isreal (  d  )  ||  d < 0  ||  1 < d
    
    % New value is no good , restore old string
    e.String = e.UserData.str ;
    
  % New value is fine
  else
    
    % Store new value
    e.UserData.str = e.String ;
    e.UserData.num = d ;
    
  end % check
  
  
end % ecutoff_cb


% Done button
function  bdone_cb (  b  ,  ~  )
  
  % makmergetool object
  h = b.Parent.UserData ;
  
  % Disable all push buttons
  h.bmerge.Enable = 'off' ;
  h.breject.Enable = 'off' ;
  h.bdone.Enable = 'off' ;
  
  % Deselect
  h.select( 0 )
  
  % Delete cluster graphics objects
  delete (  [ h.avg ; h.pca ; h.sel ; h.isi ; h.ph2 ; h.sh2 ]  )
  
  % Resume execution of merge method
  uiresume (  b.Parent  )
  
end % bdone_cb


% Cluster line callback
function  cline_cb (  l  ,  ~  )
  
  % Callback figure
  f = gcbf ;
  
  % makmergetool object
  h = f.UserData ;
  
  % Update selection
  h.select( l.UserData )
  
end % cline_cb

