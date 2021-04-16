
function  varargout = makprep (  varargin  )
% 
% defaults = makprep
% [ d , s ] = makprep (  raw  )
% [ d , s ] = makprep (  parameters  ,  raw  )
% [ d , s ] = makprep (  ...  ,  raw  ,  sid  )
% [ d , s ] = makprep (  ...  ,  raw  ,  eid  ,  sid  )
% makprep (  ...  ,  raw  ,  prep  ,  ...  )
% 
% MET Analysis Kit pre-processor. This function will convert an
% experiment's worth of data that was collected using MET into a new,
% Matlab-based format that is convenient for analysis. The output contains
% the results of automated spike sorting that must be manually verified.
% 
% With no input arguments, a struct containing default pre-processing
% parameters is returned. When at least one input argument is provided,
% then the function returns a struct with pre-processed data in d and the
% results of automated spike sorting in s. If prep is provided then no
% output arguments are produced ; d and s are saved in separate files,
% instead.
% 
% Automated spike sorting is achieved following the algorithms of Fee et
% al. (1996) and Hill et al. (2011 ; see UltraMegaSort2000). The outline is
% as follows:
% 
%   1) Align spikes to their peak
%   2) Principal component analysis , take the first N components that
%      explain a significant proportion of waveform variance. Optionally
%      apply a weighted window to each aligned waveform.
%   3) Initial waveform clustering. A large number of clusters are created,
%      far more than there are single units recorded on the electrode. This
%      allows small clusters to trace waveforms that change shape slowly
%      over time.
%   4) Compute the raw interface energy between each spike cluster
%   5) Estimate the connection strength cutoff point
%   6) Automatically merge spike clusters , starting with the pair with the
%      largest connection strength and moving downwards until the cutoff
%      point is reached. At this stage, small clusters following
%      irregularly shaped clouds of waveforms in component space will be
%      combined due to their strong connection strengths.
% 
% 
% Input arguments
% 
% raw - String naming the directory that contains all data for an
%   experiment. This will typically be the MET date directory and contain
%   several session directories, as well as one "nsp" directory that
%   contains raw electrophysiological recordings from a Blackrock
%   Microsystem's Cerebus system. The nsp directory will have one sub
%   directory for each session of data and contain .nev and .ns* files, one
%   set per trial. All session directories in raw must have the format
%   <Subject ID>.<Experiment ID>.<Session ID>[.<Tag1>][.<Tag2>][...], while
%   raw/nsp session directories must have the format
%   <Subject ID>_<Experiment ID>_<Session ID>[_...]. All session
%   directories must have the same Subject ID.
% 
% parameters - optional - A struct that has the same format as the defaults
%   struct, with certain values set as required.
% 
% sid - optional - If there are multiple sessions that need to be combined
%   then their session identifiers can be listed in sid. Alternatively,
%   this list can exclude any sessions that should be ignored. If left out,
%   then all sessions found in directory raw are combined. However,
%   sessions can only be combined if they are compatible (see note below).
%   Sessions will be combined by order of their ID, not the order that ID's
%   are given in sid.
% 
% eid - optional - If multiple sessions need to be combined from multiple
%   experiments, then eid is the experiment identifier for each session.
%   eid and sid must have the same number of elements. Experiments will be
%   combined by order of their ID, not the order that ID's are given in
%   eid.
% 
% prep - optional - String naming the destination directory where
%   pre-processed data will be stored in a .mat file, instead of returning
%   in an output argument. The output argument d is instead saved in a
%   file called prep/<Subject ID>.<Experiment ID>[.tags].mat, in a variable
%   called d ; the output argument s is saved in
%   prep/<Subject ID>.<Experiment ID>[.tags].spksort.mat, in a variable
%   called s. makprep will not write over existing files of the same name.
%   Output files will have all write permissions removed.
% 
% 
% The parameters struct will have the following fields:
% 
% .outcome - A set of characters that say which trial outcomes to allow.
%   The character mappings are:
%   
%   c - correct trials
%   f - failed trials
%   b - broken trials
%   i - ignored trials
%   a - aborted trials
% 
%   Thus, if .outcome = 'cf' then only correct and failed trials are kept.
%   This is the default.
% 
% .broken - A sub-struct that says what to do about broken trials. Since
%   poor eye-tracking or poor saccadic accuracy on part of the subject may
%   have led to the miss-classification of some trials as broken. Contains
%   fields:
% 
%   .reclassify - If non-zero then broken trials are reclassified as either
%     correct or failed. Default true.
%   .state - The break must have occurred following the onset of the named
%     task-logic state. Default 'reactim'.
%   .minlatency - The minimum duration from the onset of the named state to
%     the next MET mtarget signal, in seconds. Default 0.15.
%   .minlatkeep - If a broken trial can't be reclassified but does reach
%     the minimum latency then keep it anyway. This overrides the contents
%     of .outcome above in that a broken trial will be kept even if 'b' is
%     not in .outcome. Default true;
%   
%   Mis-classification might occur if one eye position was not within the
%   hit region of the intended target. Hence, if a trial meets the criteria
%   above, then the eye positions are examined to see if at least one eye
%   position was within a stimulus hit region. If so, then that stimulus is
%   taken to have been the target, and the trial is re-classified according
%   to the task logic. If each eye is in a different hit region then the
%   trial remains broken.
% 
% .rt - A sub-struct that says how to compute reaction times. Has fields:
%   
%   .outcome - Same format as defaults.outcome, saying what trial types to
%     compute the reaction time for. A trial with any other outcome
%     receives a default reaction time of 0. Since negative reaction times
%     can occur for trials when the subject responded too soon, it is
%     necessary to check trial outcome alongside reaction time to
%     distinguish whether 0 is the computed reaction time or not. Default
%     'cfb'.
%   .velocity - Velocity threshold. Any sample with a velocity above this
%     value is considered to be saccadic, in degrees per second. Default
%     30.
%   .accel - Acceleration threshold. Any sample with an accelleration above
%     this value is considered to be saccadic, in degrees per
%     second-squared. Default 8000.
%   .coef - Multiplicative coefficient that is applied to velocity and
%     acceleration thresholds. Default 1.0.
%   .lowpass - Cutoff frequency in Hertz of a lowpass butterworth filter
%     that is applied to the eye position recordings before computing
%     velocity and acceleration. If zero then the data is left un-filtered.
%     Default 50.
%   .order - Order of the lowpass butterworth filter. Default 5.
%   .fixdur - The minimum duration of fixation that must preceed the start
%     of the saccade, in seconds. Default 0.05.
%   .sacmin - The minimum duration prior to the last mtarget signal to
%     start searching for fixation, in seconds. Default 0.01.
%   .state - The onset of the named task-logic state is taken as time zero
%     for reaction time. If empty or the reference state can't be found
%     then the mstart time is used. Default 'reactim'.
%   .offsets - Two element vector of offsets giving the start and end of an
%     analysis window relative to the onset time of the reference state. If
%     given, the velocity and acceleration thresholds will be estimated by
%     fitting a log-normal curve to the eye data within this window.
%     Takes precedence over the .velocity and .accel manual thresholds.
%     Default [ -2 , 0 ].
%   .lognvel - Velocity threshold is taken as the point where lognvel of
%     the area under the log-normal curve is reached. This may be set to -1
%     to force the use of the .velocity set threshold. Default 0.975. 
%   .lognacc - Same as lognvel but for acceleration threshold.
%   .stim - String naming a task-logic stimulus. Any linked hit regions are
%     used to help determine which eye samples might have been pre-saccadic
%     fixations. If empty then this property is not used. Default 'fix'.
%   .verifythr - If non-zero then this is taken to be the number of seconds
%     prior to the final mtarget event in which to compare the velocity and
%     acceleration values against the measured thresholds. If the measured
%     threshold exceeds the maximum velocity or acceleration then an
%     alternative measured threshold is used in which the median plus
%     standard deviation is taken. If that still exceeds the peak value
%     then the set default value from .velocity or .accel is used. Default
%     0.25.
%   .halfmax - Logical flag. If non-zero then .verifythr is applied as
%     described above. It is then repurposed to define an analysis window
%     prior to the final mtarget event. The .logn*/.verifythr thresholds
%     are now used to search that window for significantly high velocity
%     and acceleration values. The maximum of these values are found. Half
%     the maximum observed velocity and half-max acceleration are used as
%     the final thresholds. What if the .logn*/.verifythr thresholds fail
%     to identify significant values? Then search the .verifythr window for
%     the maximum value and return half of that. Default true.
%   
%   Reaction time is defined as the time that the last saccade onto a
%   stimulus began, relative to the onset of a reference state. This is
%   computed by first averaging together binocular eye positions, then
%   filtering them with a low-pass butterworth filter. Velocity and
%   acceleration are computed from the averaged, filtered eye positions.
%   The velocity and acceleration thresholds are used to determine when the
%   eyes were fixating. Fixations may occur sacmin seconds before the final
%   mtarget signal, at latest. Another constraint is that the eyes must be
%   looking at a specified reference stimulus at the time of pre-saccadic
%   fixation, such as a fixation target. Lastly, the start of the saccade
%   must follow at least fixdur seconds of fixating eye positions that are
%   below both thresholds. If the start of a saccade cannot be found, for
%   instance if there is too much noise in the eye signal, then the time
%   of the final eye position within the hit region of the reference
%   stimulus is used ; the final mtarget time is used if no reference
%   stimulus is provided.
%   
%   If a reference state is provided, then an analysis window relative to
%   this is used to sample eye velocity and acceleration. Log-normal curves
%   are fit to each distribution and used to estimate the velocity and
%   acceleration thresholds. When no reference state is provided then the
%   manual values are used.
% 
%   Manual values are from:
% 
%     Rayner et al. 2007. Vision Research, 47(21), 2714-2726.
% 
% .photodiode - A sub-struct with information about retrieving photodiode
%   recordings. Has fields:
% 
%   .nsx - The file type suffix that contains photodiode information.
%     Default 'ns4'.
%   .eid - The electrode ID for the electrode that carried photodiode
%     information. If the photodiode was plugged into the Cerebus NSP
%     analogue input 1 then this is 129. Default 129.
%   .sync - The synchronising flip will have caused a photodiode
%     measurement that is a global peak, when there are no skipped frames.
%     This string says if it is a maximum ('max') or minimum ('min').
%     Default 'max'.
%   .flips - Says whether to match PTB flip times against either the
%     maximum ('max') or minimum ('min') of the photodiode waveform.
%     Default 'min'.
%   .leadpeak - If non-zero, then the peaks immediately preceeding those
%     identified according to .flips are used for timing. The identified
%     peaks are used to measure photometer voltage. For example, if the
%     photodiode trace is a series of maxima followed by minima and .flips
%     is 'min' then the minima are used to measure voltage and hence
%     greyscale while the maximum that preceeds each minimum is used to
%     time the flip.
%   .tailburnmax - Number of peaks to skip following end of trial mstop
%     signal before sampling between-trial mid-grey maxima. Default 8.
%   .tailburnmin - Number of peaks to skip following end of trial mstop
%     signal before sampling between-trial mid-grey minima. Default 4.
%   .useheader - Instead of estimating the inter-trial peak values from the
%     end of the trial, look at peaks before the trial starts. The
%     .tailburn* parameters will still apply, but in the opposite
%     direction. Set true to use header, defaults false to use tail.
%   .relaxedmatch - If we think that the last frame or two will be lost to
%     the photodiode recording for a known systematic reason, then make
%     this true and maknsp2ptb will match the first N peaks between the
%     recording and the PTB frame times. Default false.
%   .minprom - Minimum allowable prominence of peaks found by findpeaks, in
%     raw units. Default 1500.
%   .minwid - Minimum allowable width of peaks found by findpeaks, in
%     samples. Default 20.
%   .stdevs - The number of standard deviations above and below the median
%     post-trial peaks used to define classification boundaries. Default
%     20.
%   .skips - If non-zero then trials with frame skips are used, otherwise
%     they are discarded. Default true.
%   .skipn - The number of consecutive frames from the start and end of the
%     trial, together, with no frame skip that are required to estimate the
%     regression. Default 30.
%   .usecal - If a trial contains frame skips but .skipn consecutive frames
%     without an intervening skip cannot be found, then the mcalibrate MET
%     signals produced at the start of the trial and all following MET
%     signals can be used to estimate the regression between the data
%     acquisition clock and the MET/PTB clock. This is done if usecal is
%     non-zero. Default true.
%   .numcal - The minimum number of MET signals required to estimate
%     regression of Cerebus NSP to MET/PTB clock times. Default 5.
%   .regslope - If non-zero then the regression slope is used. If zero,
%     then a slope of 1 is assumed, because time passes at approximately
%     the same rate for a pair of computers with half-decent clocks that
%     are at the same elevation and travelling at the same speed. Default
%     false.
%   
%   The purpose of the photodiode recording is to synchronise the time
%   stamps from the data acquisition system with the timestamps recorded by
%   MET (PTB-style time stamps). This is done by locating screen flip times
%   in the photodiode recording and regressing these against the reported
%   PTB screen flip times. The line of best fit is used to convert spike
%   times from the data acquisition clock to the MET/PTB clock.
%   
%   The photodiode will have recorded from a square in the upper-left
%   corner of the stimulus screen that is mid-grey between trials, and
%   alternates between light and dark grey during a trial ; a single black
%   square is shown to indicate the initial synchronising flip in each
%   trial. To identify which maxima and minima in the recording are from
%   screen flips that occurred during the trial, it is necessary to find
%   bounds on the value of the mid-grey flips that occurred before and
%   after the trial ran. These are sampled from the tail of each trial's
%   recording. The classification boundaries are estimated by adding and
%   subtracting stdevs standard deviations above and below the median. Any
%   peak that falls outside this range will either be a light- or dark-grey
%   flip during the trial or the initial synchronising flip. The first flip
%   to exceed the mid-grey classification boundary according to sync is
%   taken to be the synchronising flip, and the following peaks are taken
%   to be in-trial flips.
%   
%   If there are frame skips then the PTB reports are used to find
%   consecutive stretches of frames at the start and end of the trial that
%   have no skips ; these are regressed against corresponding photodiode
%   peak times. When there are too few consecutive frames, the mcalibrate
%   times can optionally be used instead ; these suffer from additional
%   temporal jitter caused by the time required for MET to generate the
%   corresponding labelled TTL pulses that were recorded by the data
%   acquisition system.
% 
% .elect - A sub-struct that says how to evaluate which electrode data to
%   keep. Has fields:
%   
%   .keepall - If non-zero then all electrodes are kept and the following
%     .elect parameters are ignored. Default false.
%   .reg - String regular expression that the electrode label must match in
%     order to keep. If empty '' string then this criterion is not applied.
%     Default ''.
%   .state1 - Name of task-logic state used to align analysis window 1.
%     Uses mstart time if empty string given. Default 'present'.
%   .offset1 - Two-value vector containing the start and end times of
%     window 1 relative to the onset of the task logic state named by
%     state1, in seconds. Default [ -0.5 , 0 ].
%   .state2 , .offset2 - Same again for window 2. If state2 is empty string
%     then mstop time is used. Defaults 'present' and [ 0.05 , 0.55 ]. 
%   .rex - Regular expression string used to extract the array id and
%     channel id from the electrode label.
%     Default '^elec(?<array_id>\d+)-(?<channel_id>\d+)'
%     
%   The number of spikes observed in each window is computed for each
%   electrode on each trial. If the spike count is significantly higher in
%   window 2 than it is in window 1 then an electrode is deemed to have
%   responsive neurones. Otherwise, the electrode and its spike events are
%   rejected. The last occurence of each reference state is looked for.
% 
% .epoch - A sub-struct that defines the analysis epoch of each trial. Any
%   MET data that lies outside of this epoch is discarded. Has fields:
%   
%   .state1 - Names the task-logic state that references the start of the
%     epoch. If this is an empty string, or the state name is not in the
%     trial's task logic, then the mstart MET signal time i.e. the start of
%     the trial is used. Default 'holdfix'.
%   .offset1 - The offset in seconds from the start of state1 that the
%     epoch begins. If state1 is empty then this is ignored. Default 0.
%   .state2 - Names the task-logic state that references the end of the
%     epoch. If this is an empty string, or the state name is not in the
%     trial's task logic, then the mstop MET signal time i.e. the end of
%     the trial is used. Default ''.
%   .offset2 - The offset in seconds from the start of state2 that the
%     epoch ends. If state2 is empty then this is ignored. Default 0.
%   
%   The last occurrence of each named task-logic state in a trial is used.
%   For example, in the oddoneout.txt task logic, the holdfix state can
%   occur many times if the subject makes and breaks fixation in a rapid
%   sequence ; only the final holdfix state will be used.
% 
% .spksort - A sub-struct with information on how to do the spike sorting.
%   Has fields:
%   
%   .skip - If non-zero then automated spike sorting is not performed and
%     output argument s will be empty [] or no .spksort.mat file will be
%     written.  Default false.
% 
%   .coef_int2uv - Unit conversion coefficient. Waveforms are returned from
%     .nev files as int16. To put waveforms into micro-volts, they must
%     first be converted to floating point numbers and then multiplied by
%     this coefficient. Default 1 / 4.
%   
%   .noise_thr_auv - A noise detection threshold, in absolute micro-volts.
%     This can be used to detect and discard high-amplitude noise
%     artefacts. Otherwise, initial spike clustering can be disrupted. Any
%     spike waveform with an absolute peak that reaches or exceeds this
%     value is discarded prior to any spike sorting steps. Set to Inf to
%     keep all spikes. Default 750.
%   
%   .prethr - Number of samples taken prior to threshold crossing. Default
%     12.
%   
%   .peakwin - The duration of a window following the threshold crossing to
%     look for the spike waveform peak, in seconds. Default 2e-4.
%   
%   .fs - Sampling rate of spike waveform snippet, in Hertz. Default 30000.
%   
%   .comwin - Centre-of-mass window. A vector of sample indices relative to
%     the peak. Hence sample 0 is the peak, negative values are before the
%     peak and positive are after. For calculating the waveform's centre of
%     mass, for aligment. Default -2 : 2.
%   
%   .pcawin - Apply Gaussian shaped window to aligned waveforms prior to
%     principal component analysis. The Gaussian is centred on the aligned
%     peaks of the spikes, and it is wide enough so that a percentage of
%     the area under the curve is reached by the tail end of the spike
%     waveform window. The sum of weights is equal to the number of samples
%     in one spike waveform, to preserve the scale of waveforms. If
%     non-zero, then the weighted window is applied. Default true.
%   
%   .pcaprob - The percentage of area under the Gaussian window reached by
%     the tail end of the spike waveform window, divided by 100%. This is
%     then the fraction of the area reached, which is equivalent to the
%     probability of a normal distribution. Default 0.841, approximately 1
%     standard deviation.
%   
%   .percvar - The minimum percentage of variance to be explained by
%     principal components. The first N components are computed for each
%     waveform when it takes N components to explain at least this much of
%     the waveform variance. Default 95.
%   
%   .bisecs - The number of cluster bisections to perform during initial
%     clustering. Default 6.
%   
%   .assign - The maximum number of times to assign waveforms to clusters
%     following a bisection. Each cluster's centre is recomputed following
%     each assignment. If spikes are not reassigned before the maximum
%     number of assignments is reached then the next bisection is
%     performed. Default 5.
%   
%   .minspk - The minimum number of spikes required in each cluster.
%     Default 10.
%   
%   .defcut - The default cutoff value, in connection strength, that is
%     used to terminate automated spike cluster merging. This value is used
%     for all electrodes if it is greater than zero. If it is zero then the
%     following three parameters are used in an attempt to estimate the
%     cutoff value for each electrode, separately.
% 
%   .nboot - The number of bootstrap samples to take when estimating the
%     connection-strength cutoff value. Default 2e3.
%
%   .alpha - The bootstrap confidence interval alpha value to use when
%     estimating the connection-strength cutoff value. Default 0.01.
%   
%   .ptile - The percentile to take from each bootstrap sample. A bootstrap
%     distribution of .nboot percentiles is built up and used to estimate
%     confidence intervals at the .alpha level. Default 85.
% 
% 
% Output arguments
% 
% d - Struct with pre-processed data. In most cases, data that lies outside
%   of the analysis epoch is discarded. The smallest possible numeric types
%   are used. This assumes that trials will be under a minute long. Has
%   fields:
%   
%   
%   Session meta data -
%   
%   .preppar - A copy of the makprep parameter struct that was used
%   .subject - Subject's ID string
%   .sd - Struct vector with one session descriptor for each session
%   .header - Struct vector with header information for each session
%   .footer - Struct vector with footer information for each session
%   
%     .sd, .header, and .footer are ordered the same way as the sessions
%     are appended to each other
% 
%   .tasknames - Cell vector of strings listing the names of all tasks used
%   .logicnames - Cell vector of strings listing the names of all task
%     logics used
%   .stimlinknames - Cell vector of strings listing the unique set of all
%     stimulus link names from all tasks
%   .numtrials - The total number of trials stored in d , uint16
%   .numtrodes - The total number of electrodes , uint8
%   .electrodes - A vector listing all electrode IDs , uint8
%   .thresholds - A vector listing thresholds for each electrode , in
%     register with .electrodes , single floating point
%   .elec2chan - A uint16 vector mapping electrode ID to channel ID. Give
%     the electrode ID as the index and get the channel in return.
%   .elec2probe - A uint8 vector mapping electrode ID to probe ID.
%   
%   
%   Trial meta data -
%   
%   .experiment_id - Vector of experiment id's for each trial , uint16
%   .session_id - Vector of session id's for each trial , uint8
%   .block_id - Vector of block id's for each trial , uint16
%   .trial_id - Vector of trial id's for each trial , uint16
%   .task_ind - .tasknames index of the task used for each trial , uint8
%   .logic_ind - .logicnames index of the task logic used for each trial ,
%     uint8
%   .var - Struct with one field for each task variable. Field names are
%     taken from the task variable names. Each field contains a double
%     floating point vector containing the value of that task variable on
%     each trial.
%   .hitregion - Struct with hit region data. Each field has a cell vector
%     containing the data for each trial ; data from the start of the trial
%     to the end of the analysis epoch are kept. Has fields:
%     
%     .time - Seconds from the start of the analysis epoch that each hit
%       region update occurred
%     .stimlink - .stimlinknames index of the stimulus links represented by
%       the hit regions , column vector where rows are in register with
%       rows in .hitregion
%     .hitregion - Cell array of hit regions for each stimulus link. Rows
%       are indexed by stimulus link. Columns are indexed by hit region
%       update ; this will be empty for stimulus links that use the same
%       hit region as before. When not empty, an element of the cell array
%       contains a single floating point matrix defining a set of hit
%       regions associated with one stimulus link. The format of a hit
%       region matrix is the same as defined in metctrlconst and hence the
%       MCC variable in mak/MCC.mat. In brief, a matrix defines a set of
%       either circular or rectangular areas on the screen where objects
%       associated with the stimulus link were placed. The subject could
%       look in any one of these areas in order to select that stimulus
%       link, and hence the linked task stimulus. Each row in the matrix
%       defines a different area. Circular hit regions are defined with a
%       six-column matrix, with column order:
%         [ x , y , r , N/A , N/A , ignore ] where ( x , y ) is the
%       coordinate of the centre of the circle in degrees of visual field
%       from the middle of the screen, r is the circle's radius in visual
%       degrees, and ignore is zero when MET did not compare the area to
%       the subject's eye or mouse position i.e. the stimulus link was
%       ignored by MET. Rectangular hit regions are defined with an eight-
%       column matrix with column order:
%         [ x , y , w , h , r , N/A , N/A , ignore ] where ( x , y ) is the
%       centre of the rectangle, w and h are the width and height in visual
%       degrees, r is the counter-clockwise rotation of the rectangle, and
%       ignore has the same meaning as before.
%   
%   .start - Double floating point vector of MET/PTB time stamps marking
%     the starting time of each trial, in seconds
%   .duration - Single floating point vector of the duration of each trial,
%     in seconds
%   .epoch - Single floating point vector of the start of the analysis
%     epoch for each trial , in seconds from the start of the trial
%   .epoch_dur - Single floating point vector of the duration of the
%     analysis epoch for each trial , in seconds
%   .outcome - Char vector of the outcome of each trial. Outcome character
%     codes are:
%     
%     c - correct trials
%     f - failed trials
%     b - broken trials
%     i - ignored trials
%     a - aborted trials
%   
%   .reactime - Single floating point vector of the reaction time on each
%     trial , in seconds from the reference state ( see parameters
%     .rt.state above )
%   .reward - uint16 vector of the milliseconds of reward automatically
%     delivered at the end of each the trial
%   .rdtype - uint16 vector of the type of reward delivered automatically
%     at the end of each the trial
%   
%   
%   Cerebus NSP to MET/PTB clock alignment
%   
%   .clock.align_method - Char vector saying which alignment method was
%     used for each trial. For both types, robuts linear regression is
%     used ( see - help robustfit ). The resulting coefficients transform
%     Cerebus NSP time stamps (in seconds) to MET/PTB time stamps.
%     Character codes are:
%     
%     p - Frame onset times are measured from Cerebus photodiode recordings
%       and then regressed against the PTB stimulus onset time stamps. This
%       is the preferred method.
%     s - There was insufficient photodiode data. A backup method was used
%       in which mcalibrate signals from trial initialisation were recorded
%       by the Cerebus NSP digital input port. These times are regressed
%       against the corresponding MET/PTB signal times. This is less
%       accurate due to greater OS scheduling jitter in creating the TTL
%       pulses.
%   
%   .clock.yintercept - The y-intercept for the line of best fit on each
%     trial , single floating point vector
%   .clock.slope - The slope for the line of best fit on each trial ,
%     single floating point vector
%   
%   
%   Stimulus monitor PTB frame information
%   
%   .frame.time - Cell vector containing frame onset times for each trial.
%     Each element contains a single precision vector of frame onset times
%     from the start of the analysis epoch.
%   .frame.skip - Cell vector saying which frames that PTB classified as
%     having skipped , for each trial. An element contains a logical vector
%     that is true for skipped frames.
%   .frame.skipnum - uint16 vector of the total number of frame skips
%     within the analysis epoch for each trial
%   .frame.voltage = Cell vector containing the photometer voltage measured
%     for each frame. This should indicate the greyscale value of the
%     photodiode square on each frame. Any repeated values two or more
%     frames in a row indicates a skipped frame. Hence, there will be more
%     voltage measurements than frame onset times from PTB when there are
%     skipped frames.
%   .frame.classlims = A 2 x num trials matrix of classification limits for
%     identifying photometer voltage peaks, according to the value used in
%     p.photometer.flips. Any voltage peak that is greater than or less
%     than the classification limit while also being on the correct side of
%     zero is considered to be an in-trial peak that measures a frame.
%   
%   
%   Trial events , derived from MET signals mstate, mtarget, mreward, and
%   mrdtype
%   
%   .event.time - Cell vector containing event times for each trial. Each
%     element contains a single floating point vector with the time of each
%     trial's event , in seconds from the start of the analysis epoch.
%   data.event.type - Cell vector containing event type codes. Each element
%     contains a char vector indicating the type of each event. Character
%     codes are:
%     
%     s - Change of state. This happens any time that the task-stimulus
%       changes state during the course of the trial.
%     t - Change of target. This happens any time that the subject selected
%       a target on screen. The target is considered to be selected until
%       the next change of target event.
%     r - Reward delivery. This happens any time that a reward is delivered
%       automatically or manually.
%     y - Reward type. This happens any time that the type of reward is
%       switched.
%   
%   .event.data - Cell vector containing information about each specific
%     event , for each trial. Each element contains a uint16 vector with
%     information about each event in the trial. Information is contextual,
%     based on the type of event:
%     
%     Event type - Meaning of event.data (d is the data for one event)
%     
%              s - State identifier , see task logic .nstate{ d } to return
%                  the name of the state that was switched to. Task logic
%                  information is found in the session descriptor.
%              t - Target identifier , this is actually the stimulus
%                  identifier for the targeted task stimulus. The name of
%                  the task stimulus is found in the .nstim{ d } of the
%                  task logic.
%              r - Milliseconds of reward delivered
%              y - Type of reward. A index of 1, 2, ... etc
%   
%   
%   Eye data
%   
%   .eye.time - Cell vector of MET/PTB time stamps for each eye sample ,
%     for each trial. Each element is a single floating point vector of
%     times from the start of the analysis epoch.
%   .eye.gaze.leftx
%   .eye.gaze.lefty
%   .eye.gaze.rightx
%   .eye.gaze.righty -  Each of these is a cell vector containing gaze
%     positions for each trial. Each element contains a single floating
%     point vector of gaze positions in visual degrees from the middle of
%     the screen. leftx and lefty are the x- and y-coordinates of the left
%     eye ; rightx and righty are the x- and y-coordinates of the right
%     eye.
%   .eye.diameter.leftx
%   .eye.diameter.lefty
%   .eye.diameter.rightx
%   .eye.diameter.righty - The same as .gaze , but contains pupil diameters
%     for each trial , in pixels
%   
%   
%   Saccade data for the saccade associated with the reaction time
%   
%   .saccade.start - Single floating point vector of saccade onset times
%     for each trial , in seconds from the start of the trial
%   .saccade.duration - Single floating point vector of saccade duration
%     for each trial , in seconds
%   .saccade.angle - Single floating point vector of the saccade's angle on
%     each trial , in degrees. Angle is the counter-clockwise angle from
%     the rightward ( 1 , 0 ) basis vector in a Cartesian coordinate system
%     that is centred on the starting location of the saccade.
%   .saccade.amplitude - Single floating point vector of the saccade's
%     amplitude i.e. vector length on each trial , in visual degrees
%   .saccade.x_start
%   .saccade.y_start - The x- and y-coordinate of the start of each saccade
%     for each trial , in degrees from the middle of the screen
%   
%   
%   Mouse/touchscreen positions
%   
%   .mouse.time - Cell vector of mouse sample times for each trial. Each
%     element contains a single floating point vector of sample times from
%     the start of the analysis epoch , in seconds.
%   .mouse.x
%   .mouse.y - Cell vectors of mouse sample positions for each trial. Each
%     element contains a single floating point vector of the sample x- or
%     y-coordinates , in degrees from the middle fo the screen
%   
%   
%   Neurophysiology
%   
%   .spike.total - Double scalar value counting the total number of spikes
%     in the entire data set
%   .spike.count - E x N uint32 matrix contains the total number of spikes
%     observed on each electrode (rows) in each trail (columns) for all E
%     electrodes listed in .electrodes and all N trials
%   .spike.time - Cell vector of spike times on each trial. Each element
%     contains a single floating point vector of all spike times observed
%     in the trial from the start of the analysis epoch , in seconds.
%   .spike.electrode - Cell vector of each spike's electrode id , for each
%     trial. Each element contains a uint16 vector of electrode ids.
%   .spike.channel - Cell vector of each spike's channel id , for each
%     trial. Each element contains a uint16 vector of channel ids.
%   .spike.probe - Cell vector of each spike's probe id , for each trial.
%     Each element contains a uint8 vector of probe ids. E.g. different
%     Utah arrays are each assigned their own probe id.
%   .spike.unit - Cell vector of each spike's unit classification , for
%     each trial. Each element contains a uint8 vector of unit ids. This is
%     the classification that was assigned during the experiment by the
%     Cerebus NSP. For offline spikes sorting classifications, see s below.
%   
% 
% s - Struct with spike sorting data. Has fields:
%   
%   .wave_raw - Cell vector of raw spike waveforms on each electrode , as
%     returned by the Cerebus NSP. Each element contains a S x N int16
%     matrix of N waveforms (columns) each comprising of S samples (rows)
%     aligned to the threshold crossing , in arbitrary units.
%   .wave_aligned - Cell vector of waveforms aligned to their peak , for
%     each electrode. Each element contains a ( S - j ) x N single floating
%     point matrix of N waveforms (columns) each comprising of S - j
%     samples , in micro-volts. j is the number of samples from the
%     threshold crossing searched for the waveform peak , see parameter
%     .spksort.peakwin above.
%   .pca - Cell vector of principal components for each waveform on each
%     electrode. Each element contains a C x N single floating point matrix
%     of the C most significant principal components (rows) for each of N
%     waveforms (columns). Principal component analysis is performed using
%     singular value decomposition.
%   .spikes_per_cluster_init - Cell vector of spike counts per cluster
%     immediately after initial clustering , for each electrode. Each
%     element contains a uint32 vector.
%   .spikes_per_cluster - Cell vector of spike counts per cluster after
%     automatically merging clusters , for each electrode. Each element
%     contains a uint32 vector.
%   .d0 - Single floating point vector of scaling terms for each electrode
%   .E_init - Cell vector of raw interface energy matricies before any
%     cluster merging , for each electrode. Each element contains a N x N
%     upper triangular double floating point matrix of the energies between
%     each of N spike clusters , following initial clustering.
%   .E - Same as .E_init , but after automatically merging spike clusters
%   .J_cutoff - Double floating point vector of connection-strength cutoff
%     points for each electrode. These are used to determine when to stop
%     automatic spike cluster merging.
%   .mergers - Cell vector of cluster mergers for each electrode. Each
%     element contains a 2 x N uint8 matrix describing N mergers in the
%     order that they occurred. Each column says which two clusters where
%     merged , with row 1 containing the smaller cluster id and row 2
%     containing the bigger cluster id.
%   .clustind_init - Cell vector of cluster assignments for each spike
%     immediately following initial clustering , for each electrode. Each
%     element contains a uint8 vector of cluster ids.
%   .clustind - Same as .clustind_init , but following automatic cluster
%     merging.
%   
% 
% References:
% 
% Fee MS, Mitra PP, Kleinfeld D. J Neurosci Methods. 1996 Nov;69(2):175-88.
% Hill DN, Mehta SB, Kleinfeld D. J Neurosci. 2011 Jun 15;31(24):8699-705.
% UltraMegaSort2000, https://neurophysics.ucsd.edu/software.php
% 
% 
% NOTE: Session compatibility is in its infancy. For now, two sessions are
%   compatible if they have the same set of tags, the same set of task
%   variables, and if task variables have the same unique set of values.
% 
% NOTE: Another rudimentary approach is to assume that mreward signals for
%   correct performance will have identical signal times to the mstop
%   signal. These rewards are returned in d.reward and d.rdtype.
% 
% 
% Written by Jackson Smith - January 2018 - DPAG , University of Oxford
% 
  
  
  %%% Check input arguments %%%
  
  % Minimum/maximum number of input and output arguments
  nargoutchk (  0  ,  2  )
  narginchk  (  0  ,  5  )
  
  % Default parameters
  p = getdefaults ;
  
  % No input arguments , return default parameter struct
  if  nargin  ==  0  ,  varargout{ 1 } = p ;  return  ,  end
  
  % Set defaults
  prep = '' ;
  eid = [] ;
  sid = [] ;
  
  % Argument index
  i = 1 ;
  
  % Parameter struct is provided as first argument , guaranteed to have at
  % least one
  if  isstruct(  varargin{ i }  )
    
    % Make sure that all the required fields are present
    chkstruct (  p  ,  varargin{ i }  ,  'parameters'  ,  false  )
    
    % Point to struct
    p = varargin{ i } ;
    
    % Check next argument for raw directory name
    i = i + 1 ;
    
  end % parameter struct given
  
  % Look for raw data directory name
  if  i <= nargin  &&  ischar( varargin{ i } )  &&  ...
        isvector( varargin{ i } )  &&  exist( varargin{ i } , 'dir' )
      
    % Assign raw data directory name
    raw = varargin{ i } ;
    
    % Check next input arg for either pre-processing directory name or for
    % experiment/session id
    i = i + 1 ;
    
  % No valid input argument
  else
    
    error (  'MAK:makprep:raw'  ,  [ 'makprep:expect input arg %d ' , ...
      'to be raw data directory name' ]  ,  i  )
      
  end % raw dir
  
  % Raw data directory must have an nsp sub-directory
  if  ~ exist(  fullfile(  raw  ,  'nsp'  )  ,  'dir'  )
    
    error (  'MAK:makprep:raw'  ,  [ 'makprep:no nsp directory found ' ,...
      'in %s' ]  ,  raw  )
    
  end % no nsp dir
  
  % Look for pre-processing directory name
  if  i <= nargin  &&  ischar( varargin{ i } )  &&  ...
        isvector( varargin{ i } )  &&  exist( varargin{ i } , 'dir' )
      
    % Assign pre-processing data directory name
    prep = varargin{ i } ;
    
    % Check next input arg for experiment/session id
    i = i + 1 ;
          
  end % raw dir
  
  % If there are two more input arguments then this must be experiment id's
  if  nargin - 1 == i  &&  isnumeric( varargin{ i } )  &&  ...
        isvector( varargin{ i } )
      
    % Get experiment id's
    eid = varargin{ i } ;
    
    % Check last argument for session id's
    i = i + 1 ;
      
  end % experiment id's
  
  % The final argument must be session id's
  if  nargin == i  &&  isnumeric( varargin{ i } )  &&  ...
        isvector( varargin{ i } )
    
    % Get experiment id's
    sid = varargin{ i } ;
    
    % This should tick index one past the number of inputs
    i = i + 1 ;
    
  end % session id's
  
  % Check whether all input arguments were checked , i will be greater than
  % nargin if they were
  if  nargin  >=  i
    
    error (  'MAK:makprep:args'  ,  [ 'makprep:input argument ' , ...
      'format error' ]  )
    
  % Otherwise, make sure that eid has the same number of elements as sid
  elseif  ~ isempty( eid )  &&  numel( eid )  ~=  numel( sid )
    
    error (  'MAK:makprep:eid'  ,  [ 'makprep:eid and sid must have ' , ...
      'the same number of elements' ]  )
    
  end % final input checking
  
  
  %%% Constants %%%
  
  % data struct field names that will not have trials recursively removed.
  % First column gives parent field name, if applicable. Second column is
  % list of names to remove.
  SETDIF = ...
    {  ''  ,  [] ;
       ''  ,  { 'preppar' , 'subject_id' , 'sd' , 'header' , 'footer' , ...
                'tasknames' , 'logicnames' , 'stimlinknames' , ...
                'numtrials' , 'numtrodes' , 'electrodes' , 'thresholds' } ;
       'spike'  ,  { 'total' , 'count' }  } ; 
  
  % Load saved copies
  [ MC , MCC ] = makmetcon ;
  
  % MET signal name to signal identifier map
  MSID = MC.SIG' ;  MSID = struct (  MSID{ : }  ) ;
  
  % List of MET signals to search for
  CALSIG = [  MSID.mcalibrate  ,  MSID.mstart  ,  MSID.mstate  ,  ...
    MSID.mtarget  ,  MSID.mstop  ] ;
  
  % Outcome character map
  COUT = char ( zeros(  max( [ MC.OUT{ : , 2 } ] )  ,  1  ) ) ;
  
  for  i = 1 : size( MC.OUT , 1 )
    
    % Assign outcome character to each outcome value
    switch  MC.OUT{ i , 1 }
      
      case  'correct'  ,  COUT( MC.OUT{ i , 2 } ) = 'c' ;
      case   'failed'  ,  COUT( MC.OUT{ i , 2 } ) = 'f' ;
      case  'ignored'  ,  COUT( MC.OUT{ i , 2 } ) = 'i' ;
      case   'broken'  ,  COUT( MC.OUT{ i , 2 } ) = 'b' ;
      case  'aborted'  ,  COUT( MC.OUT{ i , 2 } ) = 'a' ;
      
    end % char assign
    
  end % outcome char map
  
  % Coefficients of lowpass butterworth filter
  BWCOEF = struct (  'b'  ,  []  ,  'a'  ,  []  ) ;
  
    % Filtering is requested if .lowpass is non-zero
    if  p.rt.lowpass

      [ BWCOEF.b , BWCOEF.a ] = butter (  p.rt.order  ,  ...
        p.rt.lowpass  /  ( MCC.SHM.EYE.SHZ / 2 )  ) ;

    end % get filter coefficients
    
  % Analysis epoch state name and offset
  EPOCHF = {  p.epoch.state1  ,  p.epoch.offset1  ;
              p.epoch.state2  ,  p.epoch.offset2  } ;
            
	% List of MET signal identifiers to map to event character codes
  MID2CH = [ MSID.mstate , MSID.mtarget , MSID.mreward , MSID.mrdtype ] ;
            
	% MET signal map to event characters
  CEVNT = char ( zeros(  max( [ MC.SIG{ : , 2 } ] ) + 1  ,  1  ) ) ;
  
  for  i = 1 : size( MC.SIG , 1 )
    
    % Assign event character to certain MET signals
    switch  MC.SIG{ i , 1 }
      
      case  'mstate'   ,  CEVNT( i ) = 's' ;
      case  'mtarget'  ,  CEVNT( i ) = 't' ;
      case  'mreward'  ,  CEVNT( i ) = 'r' ;
      case  'mrdtype'  ,  CEVNT( i ) = 'y' ;
      
    end % char assign
    
  end % event character map
  
  % Electrode probe (ELECPR) and channel index (ELECCH) map , from electrod
  % id number. Initialised empty , will be set after loading first .nev
  % file. The ELECFL map is a set of flags saying whether or not electrode
  % label matches p.elec.reg, and hence whether to keep it. ELECTH contains
  % a vector of electrode thresholds in micro-volts as single precision
  % floating point numbers.
  ELECPR = [] ;
  ELECCH = [] ;
  ELECFL = [] ;
  ELECTH = [] ;
  
  % Number of electrodes
  ENUM = 128 ;
  
  % We guess that this insertion reason value is associated only with
  % digital input events , so as to filter out serial input events
  NEV_INSERTREASON_DIGIN = 1 ;
  
  % Eye data column index map
  EYE = struct (  'LX' ,  1 ,  'LY' ,  2 ,  'RX' ,  3 ,  'RY' ,  4  ) ;
  
  % Get ready to search the raw directory for session directories by
  % building a regular expression string that will match session
  % directories created by MET.
  SESSREX = [  '(?<subject_id>\w+)'  ,  MCC.REX.SESSDIR  ] ;
  
  % Trial descriptor file format string: raw session dir / trials dir /
  % trial dir / parameter file
  TDFILE = fullfile ( '%s' , MC.SESS.TRIAL , '%d' , 'param_%d.mat' ) ;
  
  % MET signal file format string
  METSIG = fullfile ( '%s' , MC.SESS.TRIAL , '%d' , 'metsigs_%d.mat' ) ;
  
  % MET eye signal file format string
  EYEPOS = fullfile ( '%s' , MC.SESS.TRIAL , '%d' , 'eyepos_%d.mat' ) ;
  
  % MET hit region file format string
  HITREG = fullfile ( '%s' , MC.SESS.TRIAL , '%d' , 'hitregion_%d.mat' ) ;
  
  % PTB frame time stamp data file format string
  PTBTIM = fullfile ( '%s' , MC.SESS.TRIAL , '%d' , 'ptbframes_%d.mat' ) ;
  
  % Cerebus NSP event file format string
  NSPEVE = fullfile (  '%s'  ,  'trial_%d.nev'  ) ;
  
  % Photodiode recording file format string
  PHODIO = fullfile (  '%s'  ,  [ 'trial_%d.' , p.photodiode.nsx ]  ) ;
  
  % Output files are requested
  FNOUTB = '' ;
  FNOUTD = '' ;
  FNOUTS = '' ;
  if  ~ isempty (  prep  )
    
    % Base file name: prep directory/subject id.exp id[.tags]
    FNOUTB = fullfile (  prep  ,  '%s.%d%s'  ) ;
  
  end % output file names
  
  
  %%% Make sure that parallel pool is open %%%
  
  gcp ;
  
  
  %%% Find session directories %%%
  
  % Search raw data directory
  d = dir (  raw  ) ;
  
  % Match directory names against the session directory format , parse out
  % subject id , experiment id , session id , and tags
  rex = regexp (  { d.name }  ,  SESSREX  ,  'names'  ) ;
  
  % Find directory indices that do not have MET session names
  i = cellfun (  @isempty  ,  rex  ) ;
  
  % Discard those directory names
  d( i ) = [] ;
  
  % Combine into a struct vector
  rex = [  rex{ : }  ] ;
  
  % No session directories found
  if  isempty (  rex  )
    
    error (  'MAK:makprep:nosess'  ,  [ 'makprep:no session ' , ...
      'directories found in %s' ]  ,  raw  )
    
  end % no sess dirs
  
  % Turn all experiment and session id's into scalar doubles
  for  i = 1 : numel( rex )
    rex( i ).experiment_id = str2double (  rex( i ).experiment_id  ) ;
    rex( i ).session_id    = str2double (  rex( i ).session_id     ) ;
  end
  
  % Grab the subject id
  subid = rex( 1 ).subject_id ;
  
  % Check that all subject id's are the same
  if  any( ~ strcmp(  subid  ,  { rex.subject_id }  ) )
    
    error (  'MAK:makprep:subjectid'  ,  [ 'makprep:session dirs in ' , ...
      '%s for different subjects' ]  ,  raw  )
    
  end % subject id's
  
  % No session identifiers , get all exp and sess id's
  if  isempty (  sid  )
    sid = [  rex.session_id  ] ;
    eid = [  rex.experiment_id  ] ;
  end % session id's
  
  % No experiment id's given as input argument
  if  isempty (  eid  )
    
    % All session dirs must have the same experiment id
    if  any(  rex( 1 ).experiment_id  ~=  [ rex.experiment_id ]  )
    
      error (  'MAK:makprep:expid'  ,  [ 'makprep:all session dirs ' , ...
        'in %s must have the same experiment id' ]  ,  raw  )
    
    end % same exp id
    
    % Get exp id for each session id
    eid = zeros ( size(  sid  ) ) ;
    eid( : ) = rex( 1 ).experiment_id ;
    
  end % experiment id's
  
  % Make sure that sessions are ordered first by experiment id and second
  % by session id , numerical rather than string sorting
  [ ~ , i ] = sortrows (  [ eid( : ) , sid( : ) ]  ) ;
  eid = eid( i ) ;
  sid = sid( i ) ;
  
  
  %%% Output file names %%%
  
  % Output files are requested
  if  ~ isempty (  prep  )
    
    % Find first set of tags belonging to a session that will be imported
    i = [ rex.experiment_id ] == eid( 1 )  &  ...
        [ rex.session_id    ] == sid( 1 ) ;
    
    % Turn base format string into specific file names
    FNOUTB = sprintf (  FNOUTB  ,  ...
      subid  ,  min( eid )  ,  rex( i ).tags  ) ;
    
    % File name for pre-processed data
    FNOUTD = [  FNOUTB  ,  '.mat'  ] ;
    
    % File name for spike sorting data , format string
    FNOUTS = [  FNOUTB  ,  '.spksort.mat'  ] ;
    
    % Check whether either file exists already
    if  exist (  FNOUTD  ,  'file'  )  ||  exist (  FNOUTS  ,  'file'  )
      
      error (  'MAK:makprep:files'  ,  [ 'makprep: output file ' , ...
        'names already in use:\n  ' , FNOUTD , '\n  ' , FNOUTS ]  )
      
    end % files exist
    
  end % output file names
  
  
  %%% Session descriptors %%%
  
  % Store used session directory names in the order that they are appended
  sdir = cell ( size(  eid  ) ) ;
  
  % Load all session descriptors
  sd = cell ( size(  eid  ) ) ;
  
  % Header
  header = cell ( size(  eid  ) ) ;
  
  % And footer files
  footer = cell ( size(  eid  ) ) ;
  
  % For each requested session
  for  i = 1 : numel( eid )
  
    % Find session directory name
    j = [ rex.experiment_id ] == eid( i )  &  ...
        [ rex.session_id    ] == sid( i ) ;
      
    % Session directory name
    sdir{ i } = fullfile (  raw  ,  d( j ).name  ) ;
    
    % Load session descriptor
    fname = fullfile (  sdir{ i }  ,  MCC.SDFNAM  ) ;
    sd{ i } = loadvar (  fname  ,  'sd'  ) ;
    
    % Header file name
    fname = fullfile (  sdir{ i }  ,  'header.mat'  ) ;
    header{ i } = loadvar (  fname  ,  'h'  ) ;
    
    % Footer file name
    fname = fullfile (  sdir{ i }  ,  'footer.mat'  ) ;
    footer{ i } = loadvar (  fname  ,  'f'  ) ;
    
  end % load sd
  
  % Collapse into struct vectors
  sd = [  sd{ : }  ] ;
  header = [  header{ : }  ] ;
  footer = [  footer{ : }  ] ;
  
  % Check that all sessions are compatible
  for  i = 2 : numel( sd )
    
    % Check that tags are the same
    if  ~ all ( strcmp(  sd( 1 ).tags  ,  sd( i ).tags  ) )
      
      error (  'MAK:makprep:tags'  ,  [ 'makprep:%s session %d.%d ' , ...
        'has different tags from %d.%d (exp id).(sess id)' ]  ,  ...
          subid  ,  sd( 1 ).experiment_id  ,  sd( 1 ).session_id  ,  ...
            sd( i ).experiment_id  ,  sd( i ).session_id  )
      
    end % tags
    
    % Check that task variables are the same
    chkstruct (  sd( 1 ).var  ,  sd( i ).var  ,  ...
      sprintf( '%s session %d.%d sd.var' , ...
        sd( i ).experiment_id  ,  sd( i ).session_id )  ,  true  )
    
  end % check tags
  
  
  %%% Find nsp session directory names %%%
  
  ndir = cell ( size(  sdir  ) ) ;
  
  % For each requested session
  for  i = 1 : numel( eid )
    
    % Directory name search string
    fname = fullfile (  raw  ,  'nsp'  ,  ...
      sprintf( '%s_%d_%d*' , subid , eid( i ) , sid( i ) )  ) ;
    
    % Search for directory
    d = dir (  fname  ) ;
    
    % It is possible that if there are 10 or more sessions that we get
    % multiple names returned when we search for single digit sessions.
    % Look for the specific file name with regular expressions.
    j = regexp (  { d.name }  ,  sprintf( '\\w+?_\\d+_%d' , sid( i ) )  ) ;
    
    % Discard non-matches
    j = cellfun (  @isempty  ,  j  ) ;
    d( j ) = [] ;
    
    % Directory not found
    if  isempty (  d  )
      
      error (  'MAK:makprep:nspdir'  ,  [ 'makprep:%s session %d.%d ' , ...
        'has no nsp directory in %s' ]  ,  subid  ,  eid( i )  ,  ...
          sid( i )  ,  fullfile( raw , 'nsp' )  )
        
    % Too many with the same subject, exp id, and sess id
    elseif  1  <  numel( d )
      
      error (  'MAK:makprep:nspdir'  ,  [ 'makprep:%s session %d.%d ' , ...
        'multiple nsp directories found in %s' ]  ,  subid  ,  ...
          eid( i )  ,  sid( i )  ,  fullfile( raw , 'nsp' )  )
      
    end % no nsp sess dir
    
    % Store nsp session directory name
    ndir{ i } = fullfile (  raw  ,  'nsp'  ,  d.name  ) ;
    
  end % nsp session dirs
  
  
  %%% Load trials %%%
  
  % First get a fresh data structure
  [ data , spksort ] = getnewdata (  p  ,  sd  ,  ...
    sum( [  footer.trial_count  ] )  ,  header  ,  footer  ) ;
  
  % Get the difference in mu spike count on all trials for all channels.
  % Rows indexed by trial , columns by electrode id.
  spkdif = zeros (  data.numtrials  ,  ENUM  ) ;
  
  % Data structure trial index
  t = 0 ;
  
  % Discarded trial logical index vector
  disco = true (  data.numtrials  ,  1  ) ;
  
  % Step through sessions
  for  i = 1 : numel( eid )
    
    % Then step through trials in the session
    for  j = 1 : footer( i ).trial_count
      
      
      % Increment trial index
      t = t  +  1 ;
      
      % Report progress
      fprintf (  '\n%s, exp %d, sess %d, trial %d (%d of %d)'  ,...
        subid  ,  eid( i )  ,  sid( i )  ,  j  ,  t  ,  data.numtrials  )
      
      
      %-- Trial outcome --%
      
      % Load data required to determine or reclassify outcome , including
      % trial descriptor ...
      fname = sprintf ( TDFILE , sdir{ i } , j , j ) ;
      td = loadvar ( fname , 'td' ) ;
      
      % ... MET signals ...
      fname = sprintf ( METSIG , sdir{ i } , j , j ) ;
      msig = loadvar ( fname ) ;
      
      % ... eye signals ...
      fname = sprintf ( EYEPOS , sdir{ i } , j , j ) ;
      eyes = loadvar ( fname ) ;
      
      % ... and hit regions
      fname = sprintf ( HITREG , sdir{ i } , j , j ) ;
      hitreg = loadvar ( fname ) ;
      
        % Old versions of MET used 5 and 7 column arrays to define circular
        % and rectangular hit regions , made obsolete by newer 6 and 8
        % column arrays. Append a non-zero ignore column to the right side
        % of any old 5 or 7 column hit regions.
        hitreg.hitregion = cellfun (  ...
          @( h ) hitcompatible( MCC.SDEF.ptb.hitregion.ncols , h )  ,  ...
            hitreg.hitregion  ,  'UniformOutput'  ,  false  ) ;
      
      % Find mstart signal
      mstart = msig.sig  ==  MSID.mstart ;
      chknumsigs (  sum( mstart )  ,  1  ,  'mstart'  ,  sd( i )  ,  td  )
      
      % Locate mstop signal , its cargo value is the trial outcome code
      mstop = msig.sig  ==  MSID.mstop  &  msig.tim( mstart ) < msig.tim ;
      chknumsigs (  sum( mstop )  ,  1  ,  'mstop'  ,  sd( i )  ,  td  )
      
      % Get outcome character
      data.outcome( t ) = COUT (  msig.crg( mstop )  ) ;
      
      
      % Trial is broken but it is a candidate for reclassification
      if  data.outcome( t ) == 'b'  &&  p.broken.reclassify
        
        % Re-assess trial outcome
        [ data.outcome( t ) , msig , minlatflg ] = makreclasstrial (  ...
          p.broken ,  MCC ,  MSID , sd( i ) ,  td ,  msig ,  eyes.eye , ...
            hitreg  ) ;
        
        % Only use minlatflg if we are allowed to keep broken trials that
        % satisfy the latency criterion
        minlatflg = minlatflg  &&  p.broken.minlatkeep ;
        
        % Report
        if  data.outcome( t )  ~=  'b'
          fprintf (  ', broken reclassed to %s'  ,  data.outcome( t )  )
        elseif  minlatflg
          fprintf ( ', broken reached min latency' )
        end
        
      else
        
        % Not broken trial but don't allow processing of this trial to
        % continue if the outcome is not valid
        minlatflg = false ;
        
      end % reclass broken trial
      
      % We skip the trial if the outcome is not accepted for analysis
      if  all (  data.outcome( t ) ~= p.outcome  )  &&  ~ minlatflg
        fprintf (  ', reject - outcome %s'  ,  data.outcome( t )  )
        continue
      end
      
      
      %-- Align Cerebus NSP and MET/PTB clocks --%
      
      % Load PTB frame time stamp data
      fname = sprintf ( PTBTIM , sdir{ i } , j , j ) ;
      ptb = loadvar ( fname ) ;
      
      % PTB reported skipped frames and trials with frame skips are
      % rejected
      if  ~ p.photodiode.skips  &&  any (  ptb.Missed  )
        fprintf (  ', reject - skipped frames'  )
        continue
      end
      
      % Load Cerebus NSP event file
      fname = sprintf ( NSPEVE , ndir{ i } , j ) ;
      
      if  ~ exist (  fname  ,  'file'  )
        fprintf (  ', reject - NEV file missing'  )
        continue
      end
      
      nev = openNEV (  fname  ,  'nosave'  ) ;
      
        % Locate digital input events , as opposed to serial input. Take
        % the logical inverse so that we can delete serial events.
        notdev = nev.Data.SerialDigitalIO.InsertionReason  ~=  ...
               NEV_INSERTREASON_DIGIN ;
        
        % Keep only digin events
        nev.Data.SerialDigitalIO.TimeStamp( notdev )       = [ ] ;
        nev.Data.SerialDigitalIO.TimeStampSec( notdev )    = [ ] ;
        nev.Data.SerialDigitalIO.InsertionReason( notdev ) = [ ] ;
        nev.Data.SerialDigitalIO.UnparsedData( notdev )    = [ ] ;
      
      % Load photodiode recording
      fname = sprintf ( PHODIO , ndir{ i } , j ) ;
      
      if  ~ exist (  fname  ,  'file'  )
        fprintf (  ', reject - %s file missing'  ,  p.photodiode.nsx  )
        continue
      end
      
      nsx = openNSx (  fname  ) ;
      
      % Check whether the NSP recorded the analysis epoch
      nspepoch = chknspepoch (  MSID ,  EPOCHF ,  sd( i ) ,  td ,  nev  ) ;
      
      % Compute least-squares linear regression coefficients that transform
      % Cerebus NSP time stamps to MET/PTB time stamps. Note that since PTB
      % frame onset times are from the start of the trial, the n2pcoef will
      % convert NSP times to seconds from the start of the trial.
      [ n2pcoef , n2pvolts , n2pclim , n2perr ] = maknsp2ptb (  MSID  , ...
          nspepoch  ,  p.photodiode  ,  nev  ,  nsx  ,  ptb  ) ;
      
      % Error detected
      switch  n2perr
        
        % No error
        case  0  ,  fprintf (  ', align clocks by frame onset'  )
        
        % Unable to detect the start of the photodiode recording , reject
        % trial
        case  1  ,  fprintf (  ', reject - missing photodiode data'  )
                    continue
          
        % Unable to compute regression and mcalibrate estimate is not used,
        % reject trial
        case  2  ,  fprintf (  ', reject - unable to align nsp to ptb'  )
                    continue
          
        % Unable to compute regression , use mcalibrate time stamps to
        % estimate
        case  3
          
          % Report
          fprintf (  ', align clocks by MET signals'  )
          
          % Regress Cerebus NSP mcalibrate and following MET signals time
          % stamps against corresponding MET signal times ( seconds from
          % mstart i.e. start of the trial )
          n2pcoef = regmcalib (  p.photodiode  ,  MSID  ,  CALSIG  ,  ...
            msig  ,  nev  );
          
          % Too few MET signals available , discard trial
          if  isempty (  n2pcoef  )
            fprintf (  ', reject - too few'  )
            continue
          end
          
          % Store calibration method
          data.clock.align_method( t ) = 's' ;
          
          % Placeholder calibration limits
          n2pclim = [ 0 , 0 ] ;
        
      end % maknsp2ptb errors
      
      % Do we use the regression slope or assume a slope of 1? Here, the
      % parameter setting is to assume a slope of 1.
      if  ~ p.photodiode.regslope  ,  n2pcoef( 2 ) = 1 ;  end
      
      % Photometer measurements to volts
      data.frame.voltage{ t } = single (  n2pvolts  )  *  ...
        p.spksort.coef_int2uv ;
      data.frame.classlims( : , t ) = single (  n2pclim  )  *  ...
        p.spksort.coef_int2uv ;
      
      
      %-- Calculate reaction time --%
      
      % Compute reaction times for trials of a specified outcome
      if  any (  data.outcome( t ) == p.rt.outcome  )
        
        [ data.reactime( t ) , ...
          data.saccade.start( t ) , ...
          data.saccade.duration( t ) , ...
          data.saccade.angle( t ) , ...
          data.saccade.amplitude( t ) , ...
          data.saccade.x_start( t ) , ...
          data.saccade.y_start( t ) ] = makreactime (  MCC ,  MSID ,  ...
            BWCOEF ,  p.rt ,  sd( i ) ,  td ,  msig ,  eyes.eye ,  ...
              hitreg  ) ;
        
      end % rt
      
      
      %-- Determine analysis epoch start and end times --%
      
      % Get start and end times
      epoch = getepoch (  MSID  ,  EPOCHF  ,  sd( i )  ,  td  ,  ...
        msig  ,  ptb  ) ;
      
      % Error , epoch not fully contained by trial , reject
      if  ischar (  epoch  )
        
        fprintf (  ', reject - missing %s'  ,  epoch  )
        continue
        
      end % reject trial
      
      % Make times relative to start of trial
      epoch = epoch  -  ptb.trial_start ;
      
      % Epoch duration
      epdur = diff (  epoch  ) ;
      
      
      %-- Convert Cerebus NSP spike times to MET/PTB times --%
      
      % First , make sure that electrode probe and channel maps are set
      [ ELECPR , ELECCH , ELECFL , ELECTH ] = getelecprch (  p  ,  ...
        ENUM ,  ELECPR , ELECCH , ELECFL , ELECTH ,  nev  ) ;
      
      % Store them
      if  isempty (  data.elec2chan  )
        data.elec2chan = ELECCH ;
        data.elec2probe = ELECPR ;
      end
      
      % Compute spike times from mstart i.e. start of trial , in seconds
      nev.Data.Spikes.TimeStamp = n2pcoef  *  ...
     [  ones( size(  nev.Data.Spikes.TimeStamp  ) ) ;
        double( nev.Data.Spikes.TimeStamp )  /  ...
          double( nev.MetaTags.SampleRes )  ] ;
      
      % Then centre times on start of analysis epoch
      nev.Data.Spikes.TimeStamp = nev.Data.Spikes.TimeStamp  -  epoch( 1 );
      
      % Save as single precision floating point numbers
      nev.Data.Spikes.TimeStamp = single (  nev.Data.Spikes.TimeStamp  ) ;
      
      
      %-- Store trial's remaining data --%
      
      % Times: Start of trial is in MET/PTB time, start of analysis epoch
      % is in seconds from start of trial, all other times are in seconds
      % from start of analysis epoch
      
      % Hit regions , stored below
      
        % Centre hit region times on start of analysis epoch
        hitreg.time = hitreg.time  -  ptb.trial_start  -  epoch( 1 ) ;
        
        % Find hit regions that arrived after the end of the analysis
        % epoch ; these will be discarded. We keep those up to the end of
        % the analysis epoch, as even those that arrived before the epoch
        % may not have changed until after it began. Notice that since
        % times are from the beginning of the analysis epoch we must use
        % the epoch duration rather than the epoch's end time.
        k = epdur  <  hitreg.time ;
        
        % Discard hit regions arriving after the analysis epoch
        hitreg.time( k ) = [] ;
        hitreg.hitregion( k , : ) = [] ;
        
        % Determine the index of each link's name , in the
        % data.stimlinknames list. Return as column vector
        hitreg.stimlink = cellfun (  ...
          @( n ) uint8( sfind( n , data.stimlinknames ) )  ,  ...
            hitreg.stimlink  )' ;
        
        % Hit region list transposed and converted to single precision
        % floating point
        hitreg.hitregion = cellfun (  @single  ,  hitreg.hitregion'  ,  ...
          'UniformOutput'   ,  false  ) ;
      
      % Trial's meta-data
      data.experiment_id( t ) = sd( i ).experiment_id ;
      data.session_id( t ) = sd( i ).session_id ;
      data.block_id( t ) = td.block_id ;
      data.trial_id( t ) = td.trial_id ;
      data.task_ind( t ) = sfind (  data.tasknames  ,  td.task  ) ;
      data.logic_ind( t ) = sfind (  data.logicnames  ,  td.logic  ) ;
      data.var{ t } = td.var ;
      data.hitregion{ t } = hitreg ;
      data.start( t ) = ptb.trial_start ;
      data.duration( t ) = ptb.trial_end  -  ptb.trial_start ;
      data.epoch( t ) = epoch( 1 ) ;
      data.epoch_dur( t ) = epdur ;
      % data.outcome set above
      % data.reactime set above
      [ data.reward( t ) , ...
        data.rdtype( t ) ] = getreward (  MSID  ,  msig  ) ;

      % Cerebus NSP to MET/PTB clock alignment
      % data.clock.align_method set above
      data.clock.yintercept( t ) = n2pcoef( 1 ) ;
      data.clock.slope( t ) = n2pcoef( 2 ) ;

      % Stimulus monitor frame information %
        
        % Convert from microseconds since start of trial to seconds since
        % start of analysis epoch
        ptb.StimulusOnsetTime = ...
          single (  ptb.StimulusOnsetTime  )  /  1e6  -  epoch( 1 ) ;
        
        % Find frame onsets that occurred during analysis epoch. Again,
        % since times are from the start of the epoch, 0 is the time that
        % the epoch starts and the duration is the time that it ends.
        k =     0  <=  ptb.StimulusOnsetTime  &  ...
            epdur  >=  ptb.StimulusOnsetTime ;
        
        % Keep only frames times and skip flags during analysis epoch , but
        % count total number of frame skips over whole trial
        data.frame.time{ t } = ptb.StimulusOnsetTime( k )' ;
        data.frame.skip{ t } = ptb.Missed( k )' ;
        data.frame.skipnum( t ) = sum (  ptb.Missed( k )  ) ;

      % Trial events , derived from MET signals %
      
        % Centre event times on start of analysis epoch
        msig.tim = single(  msig.tim  -  ptb.trial_start  -  epoch( 1 )  );
        
        % Find signals that we want to convert into event character codes
        % that occurred during the analysis epoch. Remember, times are from
        % start of epoch so 0 and duration are start and end times.
        k = arrayfun (  @( s ) any( s  ==  MID2CH )  ,  msig.sig  ) ;
        k = k  &  0 <= msig.tim  &  msig.tim <= epdur ;
        
        % Keep event times, map to event character codes, and keep signal
        % cargo values
        data.event.time{ t } = msig.tim( k )' ;
        data.event.type{ t } = CEVNT( msig.sig(  k  )' + 1 ) ;
        data.event.data{ t } = uint16 (  msig.crg( k )'  ) ;

      % Eye data %
      
        % Centre sample times on start of analysis epoch and convert to
        % single precision floating point
        eyes.eye.time = ...
          single (  eyes.eye.time  -  ptb.trial_start  -  epoch( 1 )  ) ;
        
        % Find samples that occurred during analysis epoch. Since times are
        % from start of epoch, 0 and epoch duration are start and end
        % times.
        k = 0 <= eyes.eye.time  &  eyes.eye.time <= epdur ;
        
        % Retrieve these samples , convert to degrees/pixels , in single
        % precision format
        data.eye.time{ t } = eyes.eye.time( k )' ;
        data.eye.gaze.leftx{ t } = ...
          geteyesample (  eyes.eye.position  ,  k  ,  EYE.LX  ,  100  ) ;
        data.eye.gaze.lefty{ t } = ...
          geteyesample (  eyes.eye.position  ,  k  ,  EYE.LY  ,  100  ) ;
        data.eye.gaze.rightx{ t } = ...
          geteyesample (  eyes.eye.position  ,  k  ,  EYE.RX  ,  100  ) ;
        data.eye.gaze.righty{ t } = ...
          geteyesample (  eyes.eye.position  ,  k  ,  EYE.RY  ,  100  ) ;
        data.eye.diameter.leftx{ t } = ...
          geteyesample (  eyes.pupil.diameter  ,  k  ,  EYE.LX  ,  1  ) ;
        data.eye.diameter.lefty{ t } = ...
          geteyesample (  eyes.pupil.diameter  ,  k  ,  EYE.LY  ,  1  ) ;
        data.eye.diameter.rightx{ t } = ...
          geteyesample (  eyes.pupil.diameter  ,  k  ,  EYE.RX  ,  1  ) ;
        data.eye.diameter.righty{ t } = ...
          geteyesample (  eyes.pupil.diameter  ,  k  ,  EYE.RY  ,  1  ) ;

      % Saccade data set above

      % Mouse positions %
      
        % Centre times on start of analysis epoch and convert to single
        % floating point
        eyes.mouse.time = ...
          single (  eyes.mouse.time  -  ptb.trial_start  -  epoch( 1 )  ) ;
        
        % Find mouse samples within analysis epoch. 0 and epoch duration
        % are start and end times of epoch because sample times are now
        % from the start of the epoch.
        k = 0 <= eyes.mouse.time  &  eyes.mouse.time <= epdur ;
        
        % Get samples within analysis epoch , convert samples to singles in
        % degrees
        data.mouse.time{ t } = eyes.mouse.time( k )' ;
        data.mouse.x{ t } = ...
          geteyesample (  eyes.mouse.position  ,  k  ,  EYE.LX  ,  100  ) ;
        data.mouse.y{ t } = ...
          geteyesample (  eyes.mouse.position  ,  k  ,  EYE.LY  ,  100  ) ;

      % Spikes %
      
        % Count the difference in spikes between two analysis windows to
        % see whether electrode recorded responsive units. Note that at
        % this stage, all time stamps are relative to the start of the
        % analysis epoch, including MET signal and spike times.
        spkdif( t , : ) = getspkdiff (  p.elect , ENUM ,  MSID ,  ...
          sd( i ) ,  td ,  msig ,  nev  ) ;
        
        % Find spikes that were within analysis epoch. Remember that spike
        % times are from the start of the epoch, therefore the start and
        % end times of the epoch are 0 and the epoch duration.
        k =     0  <=  nev.Data.Spikes.TimeStamp  &  ...
            epdur  >=  nev.Data.Spikes.TimeStamp ;
          
        % We can't find the spike counts until we know which electrodes
        % we're keeping. Instead, get all spikes within the analysis
        % window, along with all information about each spike.
        data.spike.time{ t } = nev.Data.Spikes.TimeStamp( k ) ;
        data.spike.electrode{ t } = nev.Data.Spikes.Electrode( k ) ;
        data.spike.channel{ t } = ELECCH (  data.spike.electrode{ t }  ) ;
        data.spike.probe{ t } = ELECPR (  data.spike.electrode{ t }  ) ;
        data.spike.unit{ t } = int8 (  nev.Data.Spikes.Unit( k )  ) ;
        
        % Get associated waveforms , but do nothing with them until we've
        % discarded electrodes
        if  ~ p.spksort.skip
          
          spksort.wave_raw{ t } = nev.Data.Spikes.Waveform( : , k ) ;
          
        end
        
       
       %-- This is a kept trial --%
       
       disco( t ) = 0 ;
       
      
    end % trials
    
  end % sessions
  
  
  %%% Discard empty place holders %%%
  
  % For trials that were discarded, get rid of place holders in the data
  % structure
  
   % Session meta data
  data.numtrials( 1 ) = numel (  disco  )  -  sum (  disco  ) ;
  
  % Recursively discard data from structure
  data = recdisco (  data  ,  disco  ,  SETDIF  ,  2  ) ;
  
  % Spike sorting information
  if  ~ p.spksort.skip
    
    spksort.wave_raw( disco ) = [] ;
    
  end
  
  
  %%% Consolidate structs %%%
  
  % Report
  fprintf (  '\nConsolidating task variables and hit regions\n'  )
  
  % At this point, data.var and data.hitreg are cell arrays of individual
  % structs. We would rather have structs where each field contains a
  % vector. Consolidate both fields this way. hitreg's fields have cell
  % vectors, var's have numeric vectors.
  data.hitregion = consolidate (  data.hitregion  ,  false  ) ;
  data.var    = consolidate (  data.var     ,  true   ) ;
  
  
  %%% Discard high-amplitude noisy spikes %%%
  
  % Do this only if the threshold is finite
  if  ~ isinf (  p.spksort.noise_thr_auv  )
    
    % Convert threshold from micro-volts to unitless integer, comparable to
    % the raw nev waveform data
    p.spksort.noise_thr_auv = ...
      p.spksort.noise_thr_auv  ./  p.spksort.coef_int2uv ;
    
    % Check each trial
    for  t = 1 : data.numtrials
      
      % Skip trial if there are no spikes
      if  isempty( spksort.wave_raw{ t } )  ,  continue  ,  end
      
      % Find all spike waveforms with an absolute peak that reaches or
      % exceeds the noise threshold
      j = p.spksort.noise_thr_auv  <=  ...
        max (  abs(  spksort.wave_raw{ t }  )  ,  [ ]  ,  1  ) ;

      % Discard spikes
      for  F = { 'time' , 'electrode' , 'channel' , 'probe' , 'unit' }

        data.spike.( F{ 1 } ){ t }( j ) = [] ;

      end % disco

      % Discard waveforms
      spksort.wave_raw{ t }( : , j ) = [] ;

    end % trials
    
  end % discard high-amp noise
  
  
  %%% Count spikes %%%
  
  % Report
  fprintf (  'Counting spikes\n'  )
  
  % Prepare parfor variables. X is used so that the whole of data is not
  % copied for each parpool worker.
  X1 = data.spike.electrode ;
  
  % Count spikes per trial/electrode
  parfor  ip = 1 : ENUM
    
    Y( ip , : ) = cellfun(  @( eid ) uint32( sum( ip  ==  eid ) ) ,  X1  );
    
  end % spikes per trial/elec
  
  % Save parfor output
  data.spike.count = Y ;
  
  % Count all spikes per electrode
  spe = sum (  data.spike.count  ,  2  )' ;
  
  % Count all spikes in the experiment
  data.spike.total = sum ( double(  data.spike.count( : )  ) ) ;
  
  
  %%% Discard electrodes %%%
  
  % Report
  fprintf (  ...
    'Discarding unwanted electrodes or those with no responsive units\n'  )
  
  % Perform ttest on spike count differences unless all electrodes kept ...
  if  p.elect.keepall
    
    j = ones (  1  ,  ENUM  ) ;
    
  else
    
    j = ttest (  spkdif  ) ;
    
  end
  
  % ... then set any NaN value to zero and convert to logical
  j( isnan( j ) ) = 0 ;
  j = logical (  j  ) ;
  
  % Spike sorting is being skipped , so set .minspk to - Inf so that all
  % electrodes pass test
  if  p.spksort.skip  ,  p.spksort.minspk = - Inf ;  end
  
  % Find electrodes to keep. These have the correct electrode label , they
  % recorded responsive units , and they have enough spikes to cluster.
  j = ELECFL  &  j  &  p.spksort.minspk <= spe ;
  
  % Find electrodes to discard. These lack one of the properties listed
  % above.
  i = find (  ~ j  ) ;
  
  % Remove spikes from total sum
  data.spike.total = data.spike.total  -  sum (  spe( i )  ) ;
  
  % Remove spike count per trial, per electrode
  data.spike.count( i , : ) = [] ;
  
  % Store the list of kept electrode id's
  data.electrodes = uint8 ( find(  j  ) ) ;
  
  % Number of kept electrodes
  data.numtrodes( 1 ) = numel (  data.electrodes  ) ;
  
  % Electrode thresholds
  data.thresholds = ELECTH (  j  ) ;
  
  % Step through each trial
  for  t = 1 : data.numtrials
    
    % Point to electrode id of spikes
    eid = data.spike.electrode{ t } ;
    
    % Logical index for spikes
    j = false ( size(  eid  ) ) ;
    
    % Find all spikes from discarded electrodes
    for  eid = i  ,  j = j  |  eid == data.spike.electrode{ t } ;  end
    
    % Discard spikes
    for  F = { 'time' , 'electrode' , 'channel' , 'probe' , 'unit' }
      
      data.spike.( F{ 1 } ){ t }( j ) = [] ;
      
    end % disco
    
    % Discard waveforms
    if  ~ p.spksort.skip
      
      spksort.wave_raw{ t }( : , j ) = [] ;
      
    end
    
  end % trials
  
  
  %%% Perform automated spike sorting %%%
  
  if  ~ p.spksort.skip
  
  
    %-- Group waveforms by electrode --%

    % Report
    fprintf (  'Group waveforms by electrode\n'  )

    % Up to this point, spksort.wave_raw is a 1 x data.numtrials cell
    % vector with waveforms from each trial in chronological order.
    % Afterwards, spksort.wave_raw is a data.numtrodes x 1 cell vector of
    % waveforms grouped by electrode over all trials, in chronological
    % order.
    spksort.wave_raw = groupbytrode (  data  ,  spksort.wave_raw  ) ;


    %-- Align waveforms to peak --%

    % Report
    fprintf (  'Align waveforms to peak\n'  )

    % Prepare input/output variables
    X1 = p.spksort ;
    X2 = data.thresholds ;
    X3 = spksort.wave_raw ;
    Y = cell (  data.numtrodes  ,  1  ) ;

    % Convert waveforms to single floating point numbers in micro-volts and
    % align to first peak past threshold crossing , for each trial
    for  i = 1 :  data.numtrodes

      Y{ i } = makalignspks (  X1  ,  X2( i )  ,  X3{ i }  ) ;

    end % align waveforms

    % Save parfor output
    spksort.wave_aligned = Y ;
    
    % High-amplitude noise discard threshold is active. If so then it may
    % be that high-frequency noise confuses the spline interpolation that
    % is used to align the spikes to their peaks. This can result in new
    % high-amplitude spikes that weren't pruned from the original. Find
    % these and zero them.
    if  ~ isinf (  p.spksort.noise_thr_auv  )
      
      % Electrodes
      for  i = 1 : numel( spksort.wave_aligned )
        
        % Find high-amp noise
        j = p.spksort.noise_thr_auv  <=  ...
          max (  abs( spksort.wave_aligned{ i } )  ,  [ ]  ,  1  ) ;
        
        % Set to zero
        spksort.wave_aligned{ i }( : , j ) = 0 ;
        
      end % trodes
      
    end % noise thr
    
    % Find number of samples in an aligned waveform
    j = 0 ;
    for  i = 1 : numel( spksort.wave_aligned )
      if  ~ isempty ( spksort.wave_aligned{ i } )
        j = size (  spksort.wave_aligned{ i }  ,  1  ) ;
        break
      end
    end


    %-- Compute principal components --%

    % Report
    fprintf (  'Principal component analysis\n'  )

    % Compute weighted window to apply to each waveform before pca analysis
    pcawin = makpcawin (  j  ,  p.spksort.prethr  ,  ...
      p.spksort.pcawin * p.spksort.pcaprob  )' ;

    % Parfor variables
    X1 = p.spksort.percvar ;
    X2 = spksort.wave_aligned ;
    clear  Y

    % Find waveform PCA for each kept electrode
    parfor  ip = 1 : data.numtrodes

      % Number of spikes on electrode
      j = size (  X2{ ip }  ,  2  ) ;

      % Return components for each waveform on this electrode
      Y{ ip } = makpca (  X1  ,  repmat( pcawin , 1 , j )  .*  X2{ ip }  );

    end % pca

    % Save parfor output
    spksort.pca = Y' ;


    %-- Initial clustering --%

    % Report
    fprintf (  'Initial clustering\n'  )

    % Since the clustering process relies on sampling random values ,
    % reshuffle the random number generator's seed
    rng ( 'shuffle' )

    % Parfor variables
    X1 = p.spksort ;
    X2 = spksort.pca ;
    clear  Y

    % Cluster spikes on each electrode
    parfor  ip = 1 : data.numtrodes

      [ Y1{ ip } , Y2{ ip } , Y3{ ip } ] = makspkclust(  X1 ,  X2{ ip }  );

    end % initial clustering

    % Save parfor output
    spksort.spikes_per_cluster_init = Y1' ;
    spksort.clustind_init = Y2' ;
    spksort.d0 = [  Y3{ : }  ]' ;


    %-- Initialise energy matrix --%

    % Report
    fprintf (  'Initialise energy matrices\n'  )

    % Map input and output to simple variable names
    X1 = spksort.spikes_per_cluster_init ;
    X2 = spksort.clustind_init ;
    X3 = spksort.pca ;
    X4 = spksort.d0 ;
    Y = cell (  data.numtrodes  ,  1  ) ;

    % Compute matrix for each electrode
    for  i = 1 : data.numtrodes

      Y{ i } = makenergymat (  X1{ i } ,  X2{ i } ,  X3{ i } ,  X4( i )  );

    end % energy matrices

    % Save parfor output
    spksort.E_init = Y ;


    %-- Estimate connection strength cutoffs --%

    % Report
    fprintf (  'Estimate connection-strength cutoff \n'  )

    % Parfor variables
    X1 = p.spksort ;
    X2 = spksort.E_init ;
    X3 = spksort.spikes_per_cluster_init ;
    X4 = p.spksort.defcut ;
    clear  Y

    % Estimate connection-strength cutoff values , these say when to stop
    % merging clusters
    parfor  ip = 1 : data.numtrodes

      if  X4  <=  0
        Y( ip ) = makcutoff (  X1  ,  X2{ ip }  ,  X3{ ip }  ) ;
      else
        Y( ip ) = X4 ;
      end

    end % estimate cutoff

    % Save parfor output
    spksort.J_cutoff = Y' ;


    %-- Automated spike sorting --%

    % Report
    fprintf (  'Automated spike sorting\n'  )

    % Parfor variables
    X1 = spksort.spikes_per_cluster_init ;
    X2 = spksort.clustind_init ;
    X3 = spksort.E_init ;
    X4 = spksort.J_cutoff ;
    clear  Y

    % Merge similar clusters together
    parfor  ip = 1 : data.numtrodes

      [ Y1{ ip } , Y2{ ip } , Y3{ ip } , Y4{ ip } ] = ...
        makcmerge (  X1{ ip }  ,  X2{ ip }  ,  X3{ ip }  ,  X4( ip )  ) ;

    end % merge

    % Save parfor output
    spksort.spikes_per_cluster = Y1' ;
    spksort.clustind = Y2' ;
    spksort.E = Y3' ;
    spksort.mergers = Y4' ;
  
  
  end % automated spike sorting
  
  
  %%% Return / save data %%%
  
  % Output argument d requested
  if  1  <=  nargout  ,  varargout{ 1 } = data ;  end
  
  % Output argument s requested
  if  2  <=  nargout
    
    % But spike sorting skipped
    if  p.spksort.skip
      
      % Return empty
      varargout{ 2 } = [] ;
      
    % Spike sorting done
    else
      
      % Return s
      varargout{ 2 } = spksort ;
      
    end
    
  end % out arg s
  
  
  % Pre-processing directory provided
  if  ~ isempty (  prep  )
    
    % Report
    if  p.spksort.skip
      fprintf (  'Saving data to:\n  %s\n'  ,  FNOUTD  )
    else
      fprintf (  'Saving data to:\n  %s\n  %s\n'  ,  FNOUTD  ,  FNOUTS  )
    end
    
    % Rename data
    d = data ;
    s = spksort ;
    
    % Save data
    save (  FNOUTD  ,  'd'  )
    
    % Save spike sorting data if spike sorting was done
    if  ~ p.spksort.skip  ,  save (  FNOUTS  ,  's'  )  ,  end
    
    % Remove write permissions , system command
    if  p.spksort.skip
      cstr = sprintf (  'chmod a-w %s'     ,  FNOUTD  ) ;
    else
      cstr = sprintf (  'chmod a-w %s %s'  ,  FNOUTD  ,  FNOUTS  ) ;
    end
    
    % Execute write permission removal
    if  system (  cstr  ,  '-echo'  )
      
      error (  'MAK:makprep:wperms'  ,  [ 'makprep: failed to remove ' ,...
        'write permissions via system command:\n  %s']  ,  cstr  )
      
    end % remove write permissions
    
  end % save data in pre-processing directory
  
  % Report
  fprintf (  'Done\n'  )
  

end % makprep


%%% SUB-ROUTINES %%%

% Returns a struct with the default parameters
function  defaults = getdefaults
  
  % Which trial outcomes to allow
  defaults.outcome = 'cf' ;
  
  % How to handle broken trials
  defaults.broken.reclassify = true ;
  defaults.broken.state = 'reactim' ;
  defaults.broken.minlatency = 0.15 ;
  defaults.broken.minlatkeep = true ;
  
  % Reaction times
  defaults.rt.outcome = 'cfb' ;
  defaults.rt.velocity = 30 ;
  defaults.rt.accel = 8000 ;
  defaults.rt.coef = 1.0 ;
  defaults.rt.lowpass = 50 ;
  defaults.rt.order = 5 ;
  defaults.rt.fixdur = 0.05 ;
  defaults.rt.sacmin = 0.01 ;
  defaults.rt.state = 'reactim' ;
  defaults.rt.offsets = [ -2 , 0 ] ;
  defaults.rt.lognvel = 0.975 ;
  defaults.rt.lognacc = 0.99 ;
  defaults.rt.stim = 'fix' ;
  defaults.rt.verifythr = 0.25 ;
  defaults.rt.halfmax = true ;
  
  % Photodiode processing information
  defaults.photodiode.nsx = 'ns4' ;
  defaults.photodiode.eid = 129 ;
  defaults.photodiode.sync = 'max' ;
  defaults.photodiode.flips = 'min' ;
  defaults.photodiode.leadpeak = true ;
  defaults.photodiode.tailburnmax = 8 ;
  defaults.photodiode.tailburnmin = 4 ;
  defaults.photodiode.useheader = false ;
  defaults.photodiode.relaxedmatch = false ;
  defaults.photodiode.minprom = 1500 ;
  defaults.photodiode.minwid = 0.002 ;
  defaults.photodiode.stdevs = 20 ;
  defaults.photodiode.skips = true ;
  defaults.photodiode.skipn = 30 ;
  defaults.photodiode.usecal = true ;
  defaults.photodiode.numcal = 5 ;
  defaults.photodiode.regslope = false ;
  
  % Electrode acceptance criteria
  defaults.elect.keepall = false ;
  defaults.elect.reg = '' ;
  defaults.elect.state1 = 'present' ;
  defaults.elect.offset1 = [ -0.5 , 0 ] ;
  defaults.elect.state2 = 'present' ;
  defaults.elect.offset2 = [ 0.05 , 0.55 ] ;
  defaults.elect.rex = '^elec(?<probe_id>\d+)-(?<channel_id>\d+)' ;
  
  % Define analysis epoch
  defaults.epoch.state1 = 'holdfix' ;
  defaults.epoch.offset1 = 0 ;
  defaults.epoch.state2 = '' ;
  defaults.epoch.offset2 = 0 ;
  
  % Spike sorting parameters
  defaults.spksort.skip = false ;
  defaults.spksort.coef_int2uv = 1  /  4 ;
  defaults.spksort.noise_thr_auv = 750 ;
  defaults.spksort.prethr = 12 ;
  defaults.spksort.peakwin = 2e-4 ;
  defaults.spksort.fs = 30000 ;
  defaults.spksort.comwin = -2 : +2 ;
  defaults.spksort.pcawin = true ;
  defaults.spksort.pcaprob = 0.841 ;
  defaults.spksort.percvar = 95 ;
  defaults.spksort.bisecs = 6 ;
  defaults.spksort.assign = 5 ;
  defaults.spksort.minspk = 10 ;
  defaults.spksort.defcut = 0.05 ;
  defaults.spksort.nboot = 2e3 ;
  defaults.spksort.alpha = 0.01 ;
  defaults.spksort.ptile = 85 ;
  
end % getdefaults


% Initialise a new pre-processed data structure and also a spike sorting
% data structure. Many fields have transitional data structures for
% initially loading the data. But these will change form by the end of
% pre-processing.
function  [ data , spksort ] = getnewdata ( p , sd , n , h , f )
  
  
  % Session meta data
  data.preppar = p ;
  data.subject_id = sd( 1 ).subject_id ;
  data.sd = sd ;
  data.header = h ;
  data.footer = f ;
  data.tasknames = fieldnames (  sd( 1 ).task  ) ;
  data.logicnames = fieldnames (  sd( 1 ).logic  ) ;
  data.stimlinknames = [] ;
  data.numtrials = uint16 (  n  ) ;
  data.numtrodes = zeros (  1  ,  1  ,  'uint8'  ) ;
  data.electrodes = [] ;
  data.thresholds = [] ;
  data.elec2chan = [] ;
  data.elec2probe = [] ;
  
    % Get list of all unique stim link names from all tasks , start by
    % getting task structs in a struct vector, then return their contents
    % in a cell array
    d = struct2cell ( [  sd.task  ] ) ;
    
    % Every task struct has a .link field containing a struct that may have
    % different field names for each task , so expand these sub-structs
    % into a cell array. Then get the list of field names from each struct
    % , as string row vectors.
    d = cellfun (  @( d ) fieldnames( d.link )'  ,  d  ,  ...
      'UniformOutput'  ,  false  ) ;
    
    % Expand into a single row vector of strings, and return the unique
    % list of stimulus link names
    data.stimlinknames = unique ( [  d{ : }  ] ) ;
    
  
  % Trial meta data
  data.experiment_id = zeros (  1  ,  n  ,  'uint16'  ) ;
  data.session_id = zeros (  1  ,  n  ,  'uint8'  ) ;
  data.block_id = zeros (  1  ,  n  ,  'uint16'  ) ;
  data.trial_id = zeros (  1  ,  n  ,  'uint16'  ) ;
  data.task_ind = zeros (  1  ,  n  ,  'uint8'  ) ;
  data.logic_ind = zeros (  1  ,  n  ,  'uint8'  ) ;
  data.var = cell (  1  ,  n  ) ;
  data.hitregion = cell (  1  ,  n  ) ;
  data.start = zeros (  1  ,  n  ,  'double'  ) ;
  data.duration = zeros (  1  ,  n  ,  'single'  ) ;
  data.epoch = zeros (  1  ,  n  ,  'single'  ) ;
  data.epoch_dur = zeros (  1  ,  n  ,  'single'  ) ;
  data.outcome = char ( zeros(  1  ,  n  ,  'uint8'  ) ) ;
  data.reactime = zeros (  1  ,  n  ,  'single'  ) ;
  data.reward = zeros (  1  ,  n  ,  'uint16'  ) ;
  data.rdtype =  ones (  1  ,  n  ,  'uint16'  ) ;
  
  % Cerebus NSP to MET/PTB clock alignment
  data.clock.align_method = repmat (  'p'  ,  1  ,  n  ) ;
  data.clock.yintercept = zeros (  1  ,  n  ,  'single'  ) ;
  data.clock.slope = zeros (  1  ,  n  ,  'single'  ) ;
  
  % Stimulus monitor frame information
  data.frame.time = cell (  1  ,  n  ) ;
  data.frame.skip = cell (  1  ,  n  ) ;
  data.frame.skipnum = zeros (  1  ,  n  ,  'uint16'  ) ;
  data.frame.voltage = cell (  1  ,  n  ) ;
  data.frame.classlims = zeros (  2  ,  n  ,  'single'  ) ;
  
  % Trial events , derived from MET signals
  data.event.time = cell (  1  ,  n  ) ;
  data.event.type = cell (  1  ,  n  ) ;
  data.event.data = cell (  1  ,  n  ) ;
  
  % Eye data
  data.eye.time = cell (  1  ,  n  ) ;
  data.eye.gaze.leftx = cell (  1  ,  n  ) ;
  data.eye.gaze.lefty = cell (  1  ,  n  ) ;
  data.eye.gaze.rightx = cell (  1  ,  n  ) ;
  data.eye.gaze.righty = cell (  1  ,  n  ) ;
  data.eye.diameter.leftx = cell (  1  ,  n  ) ;
  data.eye.diameter.lefty = cell (  1  ,  n  ) ;
  data.eye.diameter.rightx = cell (  1  ,  n  ) ;
  data.eye.diameter.righty = cell (  1  ,  n  ) ;
  
  % Saccade data
  data.saccade.start = zeros (  1  ,  n  ,  'single'  ) ;
  data.saccade.duration  = zeros (  1  ,  n  ,  'single'  ) ;
  data.saccade.angle = zeros (  1  ,  n  ,  'single'  ) ;
  data.saccade.amplitude = zeros (  1  ,  n  ,  'single'  ) ;
  data.saccade.x_start = zeros (  1  ,  n  ,  'single'  ) ;
  data.saccade.y_start = zeros (  1  ,  n  ,  'single'  ) ;
  
  
  % Mouse positions
  data.mouse.time = cell (  1  ,  n  ) ;
  data.mouse.x = cell (  1  ,  n  ) ;
  data.mouse.y = cell (  1  ,  n  ) ;
  
  % Spikes
  data.spike.total = [] ;
  data.spike.count = [] ;
  data.spike.time = cell (  1  ,  n  ) ;
  data.spike.electrode = cell (  1  ,  n  ) ;
  data.spike.channel = cell (  1  ,  n  ) ;
  data.spike.probe = cell (  1  ,  n  ) ;
  data.spike.unit = cell (  1  ,  n  ) ;
  
  % Spike sorting information
  spksort.wave_raw = cell (  1  ,  n  ) ;
  spksort.wave_aligned = [] ;
  spksort.pca = [] ;
  spksort.spikes_per_cluster_init = [] ;
  spksort.spikes_per_cluster      = [] ;
  spksort.d0 = [] ;
  spksort.E_init = [] ;
  spksort.E      = [] ;
  spksort.J_cutoff = [] ;
  spksort.mergers = [] ;
  spksort.clustind_init =  [] ;
  spksort.clustind      =  [] ;
  
  
end % getnewdata


% Check that new struct has all the same fields as reference. Recursive
% function ends when there are no sub-structs left. Optional contents
% check.
function  chkstruct (  ref  ,  new  ,  f  ,  chkeq  )
  
  % Loop reference field names
  for  N = fieldnames (  ref  )'  ,  n = N{ 1 } ;
    
    % Check that field is present
    if  ~ isfield (  new  ,  n  )
      
      error (  'MAK:makprep:chkstruct'  ,  [ 'makprep: ' , f , ...
        ' lacks field' , n ]  )
      
    end % check field
    
    % Sub-struct required
    if  isstruct(  ref.( n )  )
      
      % Check for sub struct
      if  isstruct(  new.( n )  )
      
        % Compare sub-structs
        chkstruct (  ref.( n ) ,  new.( n ) ,  [ f , '.' , n ] ,  chkeq  )

        % No contents check valid , continue to next param
        continue

      else

        error (  'MAK:makprep:chkstruct'  ,  [ 'makprep: ' , f , ...
          ' field' , n , ' must be a sub-struct' ]  )
        
      end
      
    end % sub-struct
    
    % Do not check equality of field contents
    if  ~ chkeq  ,  continue  ,  end
    
    % String required
    if  ischar( ref.( n ) )  &&  ...
          ( ~ ischar( new.( n ) )  ||  ~ strcmp( ref.( n ) , new.( n ) )  )
        
      error (  'MAK:makprep:chkstruct'  ,  [ 'makprep: %s.%s ' , ...
        'must be %s but it is not' ]  ,  f  ,  n  ,  ref.( n )  )
      
    % Numeric set required
    elseif  isnumeric( ref.( n ) )  &&  ( ~ isnumeric( new.( n ) )  ||  ...
        ~ isempty( setdiff(  ref.( n )  ,  new.( n )  ) ) )
        
      error (  'MAK:makprep:chkstruct'  ,  [ 'makprep: %s.%s ' , ...
        'must contain values: %s' ]  ,  f  ,  n  ,  ...
          num2str( unique(  ref.( n )  ) )  )
        
    end % compare contents
    
  end % ref fields
  
end % chkstruct


% Checks for file existence then loads specified variables. If none
% specified then all data is loaded and returned.
function  varargout = loadvar (  fname  ,  varargin  )
  
  % Cannot find file
  if  ~ exist (  fname  ,  'file'  )

    error (  'MAK:makprep:loadvar'  ,  [ 'makprep:can''t find file ' , ...
      '%s' ]  ,  fname  )

  end % no sd

  % Load variables
  v = load (  fname  ,  varargin{ : }  ) ;
  
  % No variable names given , return struct
  if  nargin  ==  1
    varargout{ 1 } = v ;
    return
  end
  
  % Otherwise , load each specified variable into each requested output
  for  i = 1 : nargout
    
    varargout{ i } = v.(  varargin{ i }  ) ; %#ok
    
  end % load output args
  
end % loadvar


% Checks hit region array to make sure it has the correct number of
% columns. If it has too few, then a column of 1's is appended to the
% right-hand side of the array. If not, then the array is returned as is.
% If an unexpected number of columns is encountered then an error is
% launched.
function  h = hitcompatible ( ncols , h )
  
  % Size of array h , [ num rows , num cols ]
  s = size (  h  ) ;
  
  % Error string
  e = '' ;
  
  % Empty matrix indicates that this stimulus used an old hit region while
  % another stimulus updated its hit region dynamically during the trial.
  % This is allowed.
  if  all (  s  ==  0  )
    
    return
    
  % This equals the expected number of columns , quit now
  elseif  any (  s( 2 )  ==  ncols  )
    
    return
    
  % Too few columns , as of MET 00.04.00 a 5 column circle hit regioin was
  % in use
  elseif  s( 2 )  <  5
    
    e = 'few' ;
    
  % Too many columns
  elseif  s( 2 )  >  max (  ncols  )
    
    e = 'many' ;
  
  end % check input
  
  % Error encountered
  if  ~ isempty (  e  )
    
    mid = sprintf (  'MAK:makprep:hitregtoo%s'  ,  e  ) ;
    error(  mid  ,  'makprep: hit region array has too %s columns'  ,  e  )
    
  end % error
  
  % Append right-hand column of ones
  h = [  h  ,  ones( s( 1 ) , 1 )  ] ;
  
end % hitcompatible


% Check the number of met signals found
function  chknumsigs (  nfound  ,  nexpected  ,  msig  ,  sd  ,  td  )
  
  if  nfound  ~=  nexpected
    
    error (  [ 'MAK:makprep:' , msig ]  ,  [ 'makprep: %d %s ' ,...
      ' signals found (%d expected) in %s, exp %d, sess %d, trial %d' ],...
        nfound  ,  msig  ,  nexpected  ,  sd.subject_id  ,  ...
        sd.experiment_id  ,  sd.session_id  ,  td.trial_id  )
  
  end
  
end % chknumsigs


% Convert nsp digital input to MET signal identifiers and cargos
function  nsp = getnspsig (  nev  )
  
  % Point to digital input data
  dio = nev.Data.SerialDigitalIO ;
  
  % Find all TTL pulses that are MET signal identifiers i.e. less than 256.
  % These must be immediately followed by a valid cargo value i.e. greater
  % than 255. Remember, cargos are uint8 bit-shifted upwards by 8 bits.
  % Signal id's are plain uint8's.
  i = find (  dio.UnparsedData( 1 : end - 1 ) < 256  &  ...
              dio.UnparsedData( 2 : end     ) > 255  ) ;
  
  % MET signals come in consecutive pairs of 16-bit TTL pulses. Take the
  % time and signal identifier from the first of each pair
  nsp.tim = dio.TimeStampSec( i ) ;
  nsp.sig = dio.UnparsedData( i ) ;
  
  % Take the cargo from the second of each pair
  nsp.crg = bitshift (  dio.UnparsedData( i + 1 )  ,  - 8  ) ;
  
end % getnspsig


% Check that NSP recorded analysis epoch , based on 16-bit TTL on the
% digital input port
function  nspepoch = chknspepoch(  MSID  ,  EPOCHF  ,  sd  ,  td  ,  nev  )
  
  % Default return value is false , epoch was not recorded
  nspepoch = false ;
  
  % Convert nsp digital input to MET signal identifiers and cargos
  nsp = getnspsig (  nev  ) ;
  
  % Get default MET signal identifiers when no state is specified
  DEFSID = [  MSID.mstart  ,  MSID.mstop  ] ;
  
  % Look for the start and end of the analysis epoch
  for  i = 1 : 2
    
    % Get parameter fields for state name and offset
    [ state , offset ] = EPOCHF{ i , : } ;
    
    % Find MET signal for this state. No state specified , look for default
    % signal.
    j = [] ;
    
    if  isempty (  state  )
      
      j = find (  nsp.sig  ==  DEFSID( i )  ,  1  ,  'first'  ) ;
     
    % State is specified and in task logic , find the last occurrence
    elseif  any ( strcmp(  state  ,  sd.logic.( td.logic ).nstate  ) )
      
      % Get state identifier , this will be the signal's cargo
      istate = sd.logic.( td.logic ).istate.( state ) ;
      
      % Find 
      j = find (  nsp.sig == MSID.mstate  &  nsp.crg == istate  ,  ...
        1  ,  'last'  ) ;
      
    end % find met signal
      
    % Can't find MET signal , return false
    if  isempty (  j  )  ,  return  ,  end
    
    % Get signal time from NSP clock
    t = nsp.tim( j ) ;
    
    % Add offset and see whether this occurs before the file starts
    % recording
    if  t  +  offset  <=  0  ,  return  ,  end
    
  end % check start and end
  
  % Epoch was recorded
  nspepoch = true ;
  
end % chknspepoch


% Robust linear regression of Cerebus NSP mcalibrate + MET signal times
% against corresponding MET signal times. Returns empty if there are too
% few signals.
function  n2pcoef = regmcalib(  par  ,  MSID  ,  CALSIG  ,  msig  ,  nev  )
  
  % Default return if there are too few calibrating signals
  n2pcoef = [] ;
  
  % Find mstart time from MET records and subtract from MET signal times
  i = msig.sig  ==  MSID.mstart ;
  msig.tim = msig.tim  -  msig.tim( i ) ;
  
  % Get NSP data
  nsp = getnspsig (  nev  ) ;
  
  % Allocate cell vectors that accumulate signal indices , by signal type
  imet = cell ( size(  CALSIG  ) ) ;
  insp = cell ( size(  CALSIG  ) ) ;
  
  % Check for each type of signal
  for  s = 1 : numel( CALSIG )
    
    % Find signals in MET and NSP records
    i = find (   nsp.sig  ==  CALSIG( s )  )' ;
    j = find (  msig.sig  ==  CALSIG( s )  )' ;
    
    % No signals , check next type
    if  isempty (  i  )  ||  isempty (  j  )  ,  continue  ,  end
    
    % mcalibrate signals will be matched by cargo value
    if  CALSIG( s )  ==  MSID.mcalibrate
      
      [ insp{ s } , imet{ s } ] = arrayfun (  ...
        @( i ) mcalmatch( nsp , msig , i , j )  ,  i  ) ;
      
    % Other signals will be matched if the same number are present in both
    % records
    elseif  numel( i )  ==  numel( j )
      
      insp{ s } = i ;
      imet{ s } = j ;
      
    end % store indices
    
  end % signal types
  
  % Collapse cells into numeric index column vectors
  insp = [  insp{ : }  ]' ;
  imet = [  imet{ : }  ]' ;
  
  % Too few signals , return empty
  if  numel( insp )  <  par.numcal  ,  return  ,  end
  
  % Compute linear regression between signal time stamps
  n2pcoef = robustfit (  nsp.tim( insp )'  ,  msig.tim( imet )  )' ;
  
end % regmcalib


% Match NSP record of mcalibrate signal to MET record. If nothing found
% then return 0's.
function  [ insp , imet ] = mcalmatch ( nsp , msig , i , j )
  
  % Find MET record that matches NSP cargo value
  k = j(  nsp.crg( i )  ==  msig.crg( j )  ) ;
  
  % No match
  if  isempty ( k )  ,  insp = 0 ;  imet = 0 ;  return  ,  end
  
  % Match
  insp = i ;
  imet = k ;
  
end % mcalmatch


% Find analysis epoch of trial , returns string naming missing state or
% state + offset on error
function  epoch = getepoch (  MSID  ,  EPOCHF  ,  sd  ,  td  ,  ...
  msig  ,  ptb  )
  
  % Default epoch start and end times
  DEFTIM = [  ptb.trial_start  ,  ptb.trial_end  ] ;

  % Epoch time vector
  epoch = [ 0 , 0 ] ;
  
  % Task logic
  l = sd.logic.( td.logic ) ;

  % Search for beginning and end of analysis epoch
  for  i = 1 : 2

    % Get parameter fields for state name and offset
    [ state , offset ] = EPOCHF{ i , : } ;

    % If task logic state name is empty, or not part of trial's task
    % logic ...
    if  isempty (  state  )  ||  ...
      ~ any( strcmp(  state ,  l.nstate  ) )

      % ... then use the default time.
      epoch( i ) = DEFTIM( i ) ;

      % Look for next time point
      continue

    end % no state found
    
    % Look for last onset of named state
    j = find (  ...
      msig.sig == MSID.mstate  &  msig.crg == l.istate.( state )  ,  ...
        1  ,  'last'  ) ;
    
    % State not found
    if  isempty (  j  )
      
      % Return state name
      epoch = state ;
      
      % And quit
      return
      
    end % state not found
    
    % Get onset time plus offset
    epoch( i ) = msig.tim( j )  +  offset ;
    
    % Check that this fits within the span of the trial
    if  (  i == 1  &&  epoch( i ) < DEFTIM( i )  )  ||  ...
        (  i == 2  &&  epoch( i ) > DEFTIM( i )  )
      
      % Return error string
      if  offset  <  0
        epoch = sprintf (  '%s - %0.4f'  ,  state  ,  offset  ) ;
      else
        epoch = sprintf (  '%s + %0.4f'  ,  state  ,  offset  ) ;
      end
      
      % Quit
      return
      
    end % epoch clips trial edge

  end % start / end epoch
  
end % getepoch


% Returns vectors that map electrode id to probe or channel numbers. Also
% test whether electrode label matches given regular expression, and
% calculates electrode threshold levels.
function  [ ELECPR , ELECCH , ELECFL , ELECTH ] = getelecprch (  p  ,  ...
  ENUM  ,  ELECPR , ELECCH , ELECFL , ELECTH  ,  nev  )
  
  % Already done , quit now
  if  ~ isempty (  ELECPR  )  ,  return  ,  end
  
  % Get electrode parameters
  elect = p.elect ;
 
  % Get spike sorting parameters
  spksort = p.spksort ;
  
  % Get list of electrode id's
  eid = [  nev.ElectrodesInfo( 1 : ENUM ).ElectrodeID  ] ;
  
  % Allocate memory to mapping vectors
  ELECPR = zeros (  size(  eid  )  ,  'uint8'   ) ;
  ELECCH = zeros (  size(  eid  )  ,  'uint16'  ) ;
  
  % Retrieve channel labels , these contain probe and channel numbers
  S = {  nev.ElectrodesInfo( 1 : ENUM ).ElectrodeLabel  } ;
  
  % Make certain that strings are row vectors
  S = cellfun (  @( s ) reshape( s , 1 , numel( s ) )  ,  S  ,  ...
    'UniformOutput'  ,  false  ) ;
  
  % Parse out probe and channel id
  r = regexp (  S  ,  elect.rex  ,  'names'  ) ;
  
  % Assign mapped values
  for  i = 1 : numel ( r )
    
    % Don't bother if no values found
    if  isempty (  r{ i }  )  ,  continue  ,  end
    
    % Convert string to numeric
    ELECPR( i ) = str2double (  r{ i }.probe_id    ) ;
    ELECCH( i ) = str2double (  r{ i }.channel_id  ) ;
    
  end % assign
  
  % See whether electrode labels match regular expression filter
  if  isempty (  elect.reg  )
    
    % Empty string means that all electrodes may be considered
    ELECFL = true ( size(  ELECPR  ) ) ;
    
  % Non-empty regular expression string
  else
    
    % Apply regular expression to labels
    r = regexp (  S  ,  elect.reg  ) ;
    
    % Raise flags on those that match
    ELECFL = ~ cellfun (  @isempty  ,  r  ) ;
    
  end % reg exp filter
  
  % Allocate threshold map
  ELECTH = zeros (  size( eid )  ,  'single'  ) ;
  
  % Calculate high threshold
  hith = single( [  nev.ElectrodesInfo( 1 : ENUM ).HighThreshold  ] )  *...
    spksort.coef_int2uv ;
  
  % Calculate low threshold
  loth = single( [  nev.ElectrodesInfo( 1 : ENUM ).LowThreshold   ] )  *...
    spksort.coef_int2uv ;
  
  % Get high thresholds
  i = 0  <  hith ;
  ELECTH( i ) = hith( i ) ;
  
  % Get low thresholds
  i = loth  <  0 ;
  ELECTH( i ) = loth( i ) ;
  
  
end % getelecprch


% Return index of string s in cell vector of strings C
function  i = sfind (  s  ,  C  )
  
  i = find ( strcmp(  s  ,  C  ) ) ;
  
end % sfind


% Find reward delivered at same time as mstop signal , this is taken to be
% the reward attached with an end state
function  [ reward , rdtype ] = getreward (  MSID  ,  msig  )
  
  % Default return values
  reward = 0 ;
  rdtype = 1 ;

  % Find mstop signal
  i = msig.sig  ==  MSID.mstop ;
  
  % mstop signal time
  mstop = msig.tim( i ) ;
  
  % Reward delivered at same time
  i = msig.sig == MSID.mreward  &  msig.tim == mstop ;
  
  % Found reward signal , return duration of reward
  if  any (  i  )  ,  reward =  msig.crg( i ) ;  end
  
  % Look for reward type at time of mstop signal
  i = msig.sig == MSID.mrdtype  &  msig.tim == mstop ;
  
  % Found one , return type of reward
  if  any (  i  )  ,  rdtype =  msig.crg( i ) ;  end
  
end % getreward


% Convert numeric format to single precision floating point then divide by
% div. Row and column indices specify what data to convert. Returns as a
% row vector.
function  s = geteyesample (  data  ,  row  ,  col  ,  div  )
  
  % Convert numeric type to double for the division , and divide to change
  % unit (centi- to full unit)
  s = double (  data( row , col )  )  /  div ;
  
  % Return single precision
  s = single (  s  ) ;
  
  % Return row vector
  s = reshape (  s  ,  1  ,  numel( s )  ) ;
  
end % geteyesample


% Count difference in number of spikes recorded on each electrode when
% comparing two analysis windows. If trial's task logic does not have
% reference states then NaN is returned. NaN should be ignored by ttest,
% unless all data is NaN then NaN is returned. This must be tested for.
function  spkdif = ...
              getspkdiff (  p ,  ENUM ,  MSID ,  sd ,  td ,  msig ,  nev  )
  
  % Default return value
  spkdif = nan (  1  ,  ENUM  ) ;
  
  % If all electrodes kept then quit now
  if  p.keepall  ,  return  ,  end
  
  % Window start and end times
  wtimes = zeros (  2  ) ;
  
  % Trial's task logic
  l = sd.logic.( td.logic ) ;
  
  % Find MET signals that are used as window references
  for  C = {  { 1 , 'state1' , 'offset1' , 'mstart' }  ,  ...
              { 2 , 'state2' , 'offset2' , 'mstop'  }  }
    
    % Get window index , state param field name , and default MET signal
    [ i , state , offset , mdef ] = C{ 1 }{ : } ;
    
    % Turn parameter field names into their contents
    state  = p.( state  ) ;
    offset = p.( offset ) ;
    
    % Empty state name means to use default MET sig
    if  isempty (  state  )
    
      j = msig.sig  ==  MSID.( mdef ) ;

    % Make sure that reference state is in trial's task logic
    elseif  any ( strcmp(  state  ,  l.nstate  ) )
      
      j = msig.sig == MSID.mstate  &  msig.crg == l.istate.( state ) ;

    % Task logic does not have reference state , return default
    else

      return
      
    end % find state
    
    % Get MET signal index , search for last occurrence of state
    j = find (  j  ,  1  ,  'last'  ) ;
    
    % Window start and end times
    wtimes( i , : ) = msig.tim( j )  +  offset ;
    
  end % find MET sig indices
  
  % Spike times
  spk = nev.Data.Spikes.TimeStamp ;
  
  % Electrode ID's
  eid = nev.Data.Spikes.Electrode ;
  
  % Find spikes in each window
  w = [ wtimes( 1 , 1 )  <=  spk  &  spk  <=  wtimes( 1 , 2 ) ;
        wtimes( 2 , 1 )  <=  spk  &  spk  <=  wtimes( 2 , 2 ) ] ;
  
  % Spike counts for each window , allocate and initialise to zero
  c = zeros (  2  ,  ENUM  ) ;
  
  % Count spikes in each window , for each electrode
  for  i = 1 : ENUM
    
    % Search for spikes from electrode i
    j = eid  ==  i ;
    
    % Count spikes in each window
    c( 1 , i ) = sum (  w( 1 , j )  ) ;
    c( 2 , i ) = sum (  w( 2 , j )  ) ;
    
  end % cound spks per eid
  
  % Return difference in spike counts
  spkdif( : ) = diff (  c  ,  1  ,  1  ) ;
  
end % getspkdiff


% Recursively discard data from structure
function  s = recdisco (  s  ,  disco  ,  sdif  ,  i  )
  
  % This is not a structure or it is empty , return now
  if  isempty (  s  )  ||  ~ isstruct (  s  )  ,  return  ,  end
  
  % Loop field names , except for a base set
  for  F = setdiff (  fieldnames( s )'  ,  sdif{ i , 2 }  ) ,  f = F{ 1 } ;
    
    % This field contains a struct , apply recursive function and continue
    % to next field , and use appropriate list of field names to remove
    if  isstruct (  s.( f )  )
      
      i = strcmp (  sdif( : , 1 )  ,  f  ) ;
      if  ~ any (  i  )  ,  i = 1 ;  end
      
      s.( f ) = recdisco (  s.( f )  ,  disco  ,  sdif  ,  i  ) ;
      
    % Vector has one element per trial , nothing to discard unused trials
    elseif  size (  s.( f )  ,  2  )  ==  numel (  disco  )
      
      s.( f )( : , disco ) = [] ;
      
    end % recursion
    
  end % field names
  
end % recdisco


% Consolidate a cell array of structs into a single struct where each field
% has a vector taken from 
function  r = consolidate (  C  ,  expand   )
  
  % First turn cell array of structs into a struct vector
  s = [  C{ : }  ]' ;
  
  % Field names
  F = fieldnames (  s  ) ;
  
  % Contents of each field , partitioned into a cell vector per field
  D = num2cell (  struct2cell( s )  ,  2  ) ;
  
  % Field by field ...
  for  i = 1 : numel( F )  ,  f = F{ i } ;
    
    % ... join field names and data vectors into a new return struct
    r.( f ) = D{ i } ;
    
    % Expand cell vector into numeric vector
    if  expand  ,  r.( f ) = [  r.( f ){ : }  ] ;  end
    
  end % new struct
  
end % consolidate


% Group waveforms by electrode rather than by trial. we is empty cell if no
% spikes recorded.
function  we = groupbytrode (  d  ,  wt  )
  
  % No spikes recorded
  if  d.spike.total  ==  0
    we = {} ;
    return
  end
  
  % Prepare parfor variables
  electrodes = d.electrodes ;
  electrode  = d.spike.electrode ;
  
  % Group waveforms
  parfor  i = 1 : d.numtrodes
    
    % Electrode id
    eid = electrodes( i ) ;
    
    % Gather spikes from each trial on this electrode
    we{ i , 1 } = cellfun (  @( w , e ) w( : , e == eid )  ,  ...
      wt  ,  electrode  ,  'UniformOutput'  ,  false  ) ;
    
    % Collapse into a single matrix
    we{ i , 1 } = [  we{ i , 1 }{ : }  ] ;
    
  end % grouping
  
end % groupbytrode

