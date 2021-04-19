
MIND Analysis Kit (MAK)

This is a toolbox for the analysis of multichannel neuroscience data.

See mak/tutorial/MAK.DataAnalysisTutorial_v1.0.pdf for a brief guide to
some of the general data analysis functions.

There are a set of functions for implementing the algorithms of Fee et al.
(1996) and Hill et al. (2011). A general suite of functions is provided for
common data processing of multi-channel spike train or LFP data, for
example, as might be required by systems neuroscientists. Some convenience
functions for creating or formatting plots are also provided.

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

makddi - Computed disparity discrimination index.

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

maklinfin - Bias-corrected linear Fisher information. Uses analytically-
  derived equations of Fisher info from Kanitscheider et al. (2015). This
  quantifies the amount of information about a stimulus discrimination that
  can be obtained by an optimal linear decoder from a set of population
  responses for one stimulus value versus another.

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

makroc - Compute ROC curves and statistics for a set of spike clusters or a
  set of time bins. That is , compute ROC statistics for the marginals of
  multivariate data sets.

makrpcorr - A wrapper for the corr( ) function that packs the RHO and PVAL
  outputs into a single N by 2 matrix. Intended for use with makfun.

makskiptime - Return start time, end time, and duration of skipped frames
  as reported by Psych Toolbox.

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


Created by Jackson Smith  -  January 2018  -  DPAG , University of Oxford
Updated by JS - April 2021 - ESI (Fries Lab)

