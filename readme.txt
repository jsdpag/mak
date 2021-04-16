
MET Analysis Kit (MAK)

This is a toolbox for the offline analysis of data collected using MET.

See mak/tutorial/MAK.DataAnalysisTutorial_v1.0.pdf for a brief guide to
some of the general data analysis functions.

The first step in offline analsis is pre-processing. This will combine like
sessions, filter out unwanted trials (e.g. broken, aborted, frame skips),
use photodiode information to convert Cerebus NSP spike times into PTB
times, reject channels that do not respond to the stimulus, and reject
spikes that are outside of some epoch of interest. After reducing the data
set down to what will be analysed, the spikes on each channel are sorted
according to the algorithms of Fee et al. (1996) and Hill et al. (2011) ;
MAK additionally applies a Gaussian shaped window to peak-aligned spike
waveforms prior to the principal component analysis. Manual checking is
required to verify which spike clusters are single cells and which are
multi-unit. Pre-processed data is saved in a single data structure, ready
for analysis.


File summary:

MC.mat - A copy of the compile-time MET constants, as returned by
  met( 'const' , 1 ). Contains a struct called 'MC'.

MCC.mat - A copy of the MET controller constants as returned by
  metctrlconst. Contains a struct called 'MCC'.

NOTE: MC.mat and MCC.mat may not contain the most current values, depending
  on the version of MET. Replace MC.mat and MCC.mat with contants returned
  by the version of MET that produced the data that is under analysis.


Main pre-processing functions:

makprep - Data pre-processor. Use this function to convert raw MET and
  Cerebus data into one compact data set. Files are filtered out by outcome
  and only data within a defined analysis epoch are kept, while electrodes
  that are unresponsive or lack a matching label are discarded. The
  smallest numerical types are used that can represent the data. There are
  two stages. First, raw data is read in, trial by trial. At this stage any
  possible outcome re-classification is done and saccade parameters are
  detected. Critically, the Cerebus times are aligned to the MET/PTB times
  using frame onset times, or MET signal times. Thus, spike times are
  aligned to trial event times, and spikes occurring within the analysis
  epoch are kept. Once trials are read in and unwanted data are discarded,
  automated spike sorting is applied to each electrode using data from all
  trials. The algorithm by Fee, Mitra, and Kleinfeld (1996) is used,
  following the UltraMegaSort2000 implementation by Daniel N Hill, but with
  the addition that a weighted window is also applied to aligned spike
  waveforms. In short, and for each electrode, singular value decomposition
  is performed on the set of waveforms. The most significant components are
  used to partition the waveforms into a large set of clusters. Clusters
  are then merged together based on a non-linear distance measure that
  prioritises pairs of dense and closely neighbouring clusters. The result
  is a set of trial event data, and another set of spike sorting data. The
  spike sorting data requires manual guidance to perform the final set of
  cluster merges.
  
  NOTE: Pre-processing requires the Parallel Computing Toolbox.

makmancmerge - Manual cluster merging. Steps through each electrode so that
  the final set of cluster merges can be assessed and chosen. Average
  waveforms, principal component spaces, inster-spike-interval
  distrubutions, and connection strengths are provided to guide manual
  merging. The result is a data set with final spike cluster assignments
  given to each spike that has not been rejected.


General pre-processing functions:

makmetcon - Returns the MET compile time and MET controller constants
  contained in MC.mat and MCC.mat.

makrepair - Used to fix errors in the naming of the raw data files. Mainly
  for data obtained when Cerebus > Central > File Storage is synchronised
  with MET using the Serial port.


Event pre-processing functions:

maknsp2ptb - Takes photodiode recording from the Cerebus system, detects
  frame onset times according to the Cerebus clock, then regresses these
  against the PTB stimulus onset times. The result are linear coefficents
  that convert Cerebus NSP time stamps to MET/PTB times that are aligned
  with event times.

makreactime - Computes the reaction time for a trial by detecting the
  subject's final saccade and measuring the latency from a given trial
  event to the start of the saccade. Saccade parameters are also returned.

makreclasstrial - Poor saccadic accuracy coupled with eye-tracker noise
  may have resulted in some trials being mis-classified. Trials with
  certain outcomes can be re-examined by this function to see if they
  should have a different outcome, which affects whether or not they should
  be discarded.


Spike sorting pre-processing functions:

makalignspks - Spike waveforms are aligned to their peak. The peak is
  detected by finding the peak sample following the threshold crossing.
  Then, a centre of mass around the peak is computed ; essentially, the
  time index is weighted by the waveform, summed, and normalised. This
  offset is applied to the peak sample time to estimate the actual peak
  time. The waveform is aligned to this time using spline interpolation.
  This is an important first step prior to singular value decomposition, as
  it will reduce the number of components required to explain waveform
  variance. This takes time.

makcind - Returns sets of linear indices for accessing portions of the
  interface enregy matrix that relate to a pair of cluster indices. This is
  mainly for adjusting the energy matrix following a merger. See makcmerge.

makcmerge - Automated spike cluster merging. Initial clusters are merged
  together based on their connection strength. The next merger is done by
  selecting the pair of clusters with the highest connection strength. This
  continues until the connection strength between clusters drops below a
  cutoff level.

makconnstrength - The raw interface energy is stored because new energies
  can easily be computed following a cluster merger. However, interface
  energy must be normalised into a connection strength, which then guides
  cluster merging. This performs such a normalisation and returns the
  connection strength matrix.

makcutoff - If requested, the connection strength cutoff can be estimated
  from the data, unless a single constant value is provided. The estimate
  is taken by taking the upper BCA boostrap confidence interval on a
  certain percentile of the inter-cluster connection strength.

makenergymat - Computes the initial raw interface energy matrix for all
  initial spike clusters. This takes much time as every single pairwise
  distance is computed between all spikes.

makmergetool - Returns a makmergetool object. This produces the GUI tool
  for manual spike merging with makmancmerge.

makpca - Performs principal component analysis on a set of waveforms using
  singular value decomposition. Only the N most significant components are
  returned that capture a given percentage of waveform variance.

makpcawin - Return a normalised segment of a Gaussian for use as a weighted
  window that is applied to each waveform prior to principal component
  analysis. If the peak of the Gaussian matches the peak of the waveform,
  then more emphasis is put on the peak of the waveform than the head or
  tail. In other words, this window helps to emphasise parts of the
  waveform with a high signal-to-noise ratio, and reduces emphasis on the
  parts with a low S:N.

makspkclust - Initial spike clustering. Partitions spikes in the PCA space
  up into a large set of small clusters. The purpose of this it to tile
  each cloud of spikes with small clusters that can be rejoined based on
  their connection strength. This allows spikes from the same unit to be
  clustered together in spite of non-linear noise, such as electrode drift.
  This matters when two non-linear clouds of spikes are close to each other
  such that linear distance based measures might cause spikes from separate
  units to be clustered. However, a non-linear connection strength will
  tend to group small clusters from the same unit, due to their density and
  proximity.


General analysis:

makbalancedz - Convert input data into balanced z-scores , as described by
  Kang and Maunsell (2012) for computing corrected choice and detect
  probabilities.

makbindvar - Bin a dependent variable by the values of an independent
  variable and then get the mean and error of each bin. This is useful for
  turning messy scatter plots into clean line plots showing the central
  tendency of the data.

makcf - Short form for cellfun( ... , 'UniformOutput' , false )

makconv - FFT based convolution providing an interface that's more
  convenient than Matlab's conv function. This can apply a kernel to a set
  of signals, removing the tailing edges of the convolution according to
  whether the kernel is predictive, symmetrical, or causal.

makcrosstalk - Calculate measures of shared signal between all pairs of
  kept electrodes. High numbers of synchronous spikes with correlated
  waveform shapes is an indication of some strong shared signal between a
  pair of channels.

makddi - Computed disparity discrimination index.

makevtim - Return the time in each trial of the last instance of a
  specified event. For instance, return the time in each trial when the
  stimulus actually appeared. This is very useful for zeroing spike times
  on trial events when there is some jitter in the event timing from trial
  to trial.

makgethits - Tests monocular eye positions against a set of stimulus hit
  regions. Flags any eye positions that fall within any hit region.

makload - Load a set of experiments into a single data structure. This will
  append manual spike sorting results to the event data. It can also add
  any auxiliary data from files written by maksave.

makfname - Return analysis data file name built for a data set returned by
  makload.

makfun - This is a powerful high-level function along the lines of cellfun.
  It can apply a set of functions to a data set that is grouped according
  to a set of variables. It can be used in combination with functions like
  bootci. Common uses include, but are by no means limited to, calculating
  psychometric curves, neural tuning curves, and peri-stimulus time
  histograms of the spiking rate for different stimulus parameters or
  behavioral outcomes.

makgabor - Implements a 1-dimensional Gabor function. Suitable for fitting
  to disparity tuning curves with Matlab's lsqcurvefit function.

makgaborfit - Finds best fitting Gabor for each curve in a set using a
  non-linear least-squares search.

makimat - Return logical index matrix in which only the upper-triangular
  half contains 'true' values. The rest are false. Logical index vectors
  can be given to futher refine which elements to return e.g. pairs of
  spike clusters with specific attributes.

makjennrich - Jennrich Test for equality between two correlation matrices.

makkeepspk - Can be used to identify spikes that are unique to one channel,
  discarding any other that is classified as being a result of cross-talk.

maklinfin - Bias-corrected linear Fisher information. Uses analytically-
  derived equations of Fisher info from Kanitscheider et al. (2015). This
  quantifies the amount of information about a stimulus discrimination that
  can be obtained by an optimal linear decoder from a set of population
  responses for one stimulus value versus another.

makmedstd - Compute median of a sample and the standard deviation around
  the median , rather than the conventional mean.

makmi - Computes empirical mutual information between a sample of signal
  values and multiple sets of output samples.

makpak - Returns specific output arguments of a given function in a single
  cell array. For use with makfun.

makpspkern - Return a convolution kernel in the shape of a postsynaptic
  potential. See Thompson, Hanes, Bichot, & Schall. 1996). J Neurophysiol
  76(6): 4040-4055.

makrccg - An implementation of Wyeth Bair's r_ccg measure of spike train
  correlation. It reveals the time-scale at which correlations occur. This
  implementation computes the pair-wise r_ccg between multiple sets of
  spike trains ; when there are multiple spike clusters available , for
  instance. It can also pre-processed data for repeat computations of
  r_ccg with different sub-sets of trials, or for bootstrapping.

makrccg2 - New implementation of r_CCG optimised for speed.

makrfgauss - Simple GUI tool for user-guided receptive field (RF)
  quantification. The user provides manual guidance in fitting a 2D
  isotropic Gaussian to RF mapping data for each electrode. First, the user
  chooses which RF mapping data set(s) to include, then provides starting
  baseline, centre and width parameters. Least-squares linear regression
  fits are then generated.

makroc - Compute ROC curves and statistics for a set of spike clusters or a
  set of time bins. That is , compute ROC statistics for the marginals of
  multivariate data sets.

makrpcorr - A wrapper for the corr( ) function that packs the RHO and PVAL
  outputs into a single N by 2 matrix. Intended for use with makfun.

maksave - Saves auxiliary data for a set of experiments. This is meant to
  be the main way that analytical results are stored.

makskiptime - Return start time, end time, and duration of skipped frames
  as reported by Psych Toolbox.

makspk - Returns spike times grouped by trial and spike cluster in a 2D
  cell array.

maksttc - Computes Cutts & Engel's STTC metric of spike train correlation
  at different time scales.

maksttc_cutts - Computes Cutts & Engel's STTC metric using their O( n^2 )
  code , for validation of maksttc.

maktiedrank - Computes tied ranks for each column of an input matrix.

makwavg - Computes weighted average of numeric data.

makxcorr - Simple re-implementation of Matlab's xcorr. Uses increased
  vectorisation for improved run times.


Plotting functions:

makaxeq - Make axes objects square and tight around the data, then match
  the x- and y-axis limits.

makrastplot - Creates a raster plot with trials along the y-axis and tick
  marks along the x-axis. Useful for visual inspection of single spike
  clusters.


Written by Jackson Smith  -  January 2018  -  DPAG , University of Oxford
