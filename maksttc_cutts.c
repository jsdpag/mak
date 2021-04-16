
/*  sttc_cutts
  
  [ sttc , DT ] = maksttc_cutts ( win , dt , A , B )
  
  A MEX adaptation of the spike_time_tiling_coefficient.c function provided
  by Cutts and Eglen at https://github. com/CCutts/
  Detecting_pairwise_correlations_in_spike_trains. This uses an O( n ^ 2 )
  algorithm. The MEX wrapper will evaluate sttc between spike trains A and 
  B within analysis window win at each delta-t value from 0 to dt, rounded
  up to the nearest millisecond. This adds another level of nested looping,
  so the overal algorithm is O( n ^ 3 ). A and B must be double floating
  point vectors of spike times in chronological order, and dt must be a
  scalar double with delta-t in seconds. win is a two-element double of the
  start and end of the window in seconds. Returns a double vector. Optional
  output DT is a list of delta-t times in register with sttc, also a double
  vector.
  
  It is added to MAK mainly for validation of maksttc.
  
  Reference: 
  
    Cutts CS, Eglen SJ. 2014. Detecting Pairwise Correlations in Spike
      Trains: An Objective Comparison of Methods and Application to the
      Study of Retinal Waves. J Neurosc, 34(43):14288-14303.
  
  Adapted by Jackson Smith - March 2018 - DPAG , University of Oxford
  
*/


/*-- Include block --*/

#include    <math.h>
#include     "mex.h"
#include  "matrix.h"


/*-- Define block --*/

#define   NARGIN  4
#define  NARGOUT  2
#define   WINARG  0
#define    DTARG  1
#define     AARG  2
#define     BARG  3
#define  STTCARG  0
#define    DTARG  1


/*** Cutts' code block ***/

/* Compute proportion of spikes from one train within delta-t of another */
double  run_P ( int  N1 , int  N2 , double  dt , double *  spike_times_1 , 
          double *  spike_times_2 )
{

  /* Variables */
	int  i ;
	int  j ;
	int  Nab ;
		
  /* Initialise */
	Nab = 0 ;
	  j = 0 ;
  
  /* Check every spike in train 1 to see if there's a spike in train 2
     within dt (don't count spike pairs) don't need to search all j each
     iteration */
	for  ( i = 0 ; i <= ( N1 - 1 ) ; i++ )
  
		while  ( j < N2 )
    {	

			if  ( fabs( spike_times_1[ i ] - spike_times_2[ j ] )  <=  dt )
      {
				Nab = Nab + 1 ;	
				break ;				
			}
			else if  ( spike_times_2[ j ]  >  spike_times_1[ i ] )
      	
				break ;
			
			else
        
				j = j + 1 ;
      
		} /* while */
	
  /* Return count , this is scaled in run_sttc */
	return  Nab ;
  
} /* run_P */


/* Compute proportion of time within delta-t of spikes */
double  run_T ( int  N1v , double  dtv , double  startv , double  endv ,
          double * spike_times_1 )
{

  /* Variables */
	double  dt = dtv ;
	double  start = startv ;
	double  end = endv ;
	   int  N1 = N1v ;
	double  time_A ;
	   int  i = 0 ;
	double  diff ;
	
  /* Maximum */
	time_A = 2  *  ( double ) N1  *  dt ;

  /* if just one spike in train */
	if  ( N1 == 1 )
  {
		
	  if  (  ( spike_times_1[ 0 ] - start )  <  dt  )
      
	    time_A = time_A - start + spike_times_1[ 0 ] - dt ;
	  
	  else if  (  ( spike_times_1[ 0 ] + dt )  >  end  )
      
	   	time_A = time_A - spike_times_1[ 0 ] - dt + end ;
	  
	}
	
  /* if more than one spike in train */
	else
  {
	
	  while  ( i < ( N1 - 1 ) )
    {
			
			diff = spike_times_1[ i + 1 ] - spike_times_1[ i ] ;
				
      /* subtract overlap */
			if  ( diff < 2 * dt )
				time_A = time_A - 2 * dt + diff ;
				
			i++;
		}
				
		/* check if spikes are within dt of the start and/or end, if so just
       need to subract overlap of first and/or last spike as all within-
       train overlaps have been accounted for */
		if  ( ( spike_times_1[ 0 ] - start )  <  dt )
		  time_A = time_A - start + spike_times_1[ 0 ] - dt ;

		if  ( ( end - spike_times_1[ N1 - 1 ] )  <  dt )
			time_A = time_A - spike_times_1[ N1 - 1 ] - dt + end ;
    
  }
	
  /* Return sum , this is scaled by run_sttc */
	return  time_A ;
  
} /* run_T */
	

/* Call this from mexFunction for each delta-t */
void  run_sttc ( int * N1v , int * N2v , double * dtv , double * Time , 
        double * index , double * spike_times_1 , double * spike_times_2 )
{

  /* Variables */
	double  TA ;
	double  TB ;
	double  PA ;
	double  PB ;
	   int  N1 = *N1v ;
	   int  N2 = *N2v ;
	double  dt = *dtv ;
	double  T ;

	
	/* Empty spike train , STTC is undefined so return NaN */
	if ( N1 == 0  ||  N2 == 0 )
  {
    *index = mxGetNaN ( ) ;
    return ;
  }
	
	/* Both trains have spikes , get duration of analysis window */
	T = Time[ 1 ] - Time[ 0 ] ;
  
  /* Compute proportion of time within delta-t of spikes in each train */
	TA = run_T ( N1 , dt , Time[ 0 ] , Time[ 1 ] , spike_times_1 ) ;
	TA = TA / T ;
  
	TB = run_T ( N2 , dt , Time[ 0 ] , Time[ 1 ] , spike_times_2 ) ;
	TB = TB / T ;
  
  /* Compute proportion of spikes from one train within delta-t of the
     other */
	PA = run_P ( N1 , N2 , dt , spike_times_1 , spike_times_2 ) ;
	PA = PA / (double) N1 ;
  
	PB = run_P ( N2 , N1 , dt , spike_times_2 , spike_times_1 ) ;
	PB = PB / (double) N2 ;
  
  /* At last , compute STTC */
	*index = 0.5  *  ( PA - TB ) / ( 1 - TB * PA )  +  
           0.5  *  ( PB - TA ) / ( 1 - TA * PB ) ;

} /* run_sttc */


/*** MEX gateway function ***/

void  mexFunction (  int nlhs  ,        mxArray * plhs[ ] ,
                     int nrhs  ,  const mxArray * prhs[ ]  )
{
  
  
  /*-- Variables --*/
  
  /* Generic counter */
  unsigned int  i = 0 ;
  
  /* delta-t value */
  double  dtv = 0 ;
  
  /* Pointers to win, A, B, sttc, and DT data */
  double  * win , * A , * B , * sttc , * DT ;
  
  /* Number of delta-t values at millisecond steps , including zero */
  unsigned int  dtn = 0 ;
  
  /* Number of spikes from A and B within limits of win */
  int  Na = 0 , Nb = 0 ;
  
  
  /*-- Input check --*/
  
  /* Must be exactly 4 input args */
  if  ( nrhs  !=  NARGIN )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:nargsin"  ,  
      "maksttc_cutts: requires %d input arguments"  ,  NARGIN  ) ;
    
  /* Must be no more than 2 output args */
  else if  ( NARGOUT  <  nlhs )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:nargsout"  ,  
      "maksttc_cutts: returns at most %d output arguments"  ,  NARGOUT  ) ;
  
  /* win must be 2 element double */
  else if  (  !mxIsDouble( prhs[ WINARG ] )  ||
               mxGetNumberOfElements( prhs[ WINARG ] ) != 2  )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:win"  ,  
      "maksttc_cutts: win must be a two-element double"  ) ;
  
  /* dt must be a scalar double */
  else if  (  !mxIsDouble( prhs[ DTARG ] )  ||
              !mxIsScalar( prhs[ DTARG ] )  )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:dt"  ,  
      "maksttc_cutts: dt must be a scalar double"  ) ;
  
  /* A and B must be doubles */
  else if  (  !mxIsDouble( prhs[ AARG ] )  ||  
              !mxIsDouble( prhs[ BARG ] )  )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:AB"  ,  
      "maksttc_cutts: A and B must be doubles"  ) ;
  
  /* Access delta-t value */
  dtv = mxGetScalar (  prhs[ DTARG ]  ) ;
  
  /* dt must not be less than zero */
  if  ( dtv  <  0 )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:neg_dt"  ,  
      "maksttc_cutts: dt must be zero or more"  ) ;
  
  
  /*-- Preparation --*/
  
  /* Number of millisecond steps from 0 to dt , rounded up */
  dtn = ceil ( dtv  /  0.001 )  +  1 ;
  
  /* Get number of spikes in total from each train */
  Na = mxGetNumberOfElements(  prhs[ AARG ]  ) ;
  Nb = mxGetNumberOfElements(  prhs[ BARG ]  ) ;
  
  /* Point to win, A, and B data */
  win = mxGetPr (  prhs[ WINARG ]  ) ;
    A = mxGetPr (  prhs[   AARG ]  ) ;
    B = mxGetPr (  prhs[   BARG ]  ) ;
    
  /* Allocate output data */
  sttc = mxCalloc (  dtn  ,  sizeof( double )  ) ;
  
  if  ( 1  <  nlhs )
    DT = mxCalloc (  dtn  ,  sizeof( double )  ) ;
    
  /* Make sure that win( 2 ) is greater than win( 1 ) */
  if  ( win[ 1 ]  <=  win[ 0 ] )
    
    mexErrMsgIdAndTxt (  "MAK:maksttc_cutts:winlim"  ,  
      "maksttc_cutts: win( 2 ) must be greater than win( 1 )"  ) ;
    
  /* Find the first spike from A that is within the window */
  while  (  *A < win[ 0 ]  &&  Na  )
  {
    A++ ;  Na-- ;
  }
    
  /* Same again for B */
  while  (  *B < win[ 0 ]  &&  Nb  )
  {
    B++ ;  Nb-- ;
  }
    
  /* Find the last spike from A in the window */
  for  (  i = 0  ;  i < Na  &&  A[ i ] <= win[ 1 ]  ;  i++  )  ;
    
  /* Subtract spikes that fall off the tail end of window */
  Na -= ( Na  -  i ) ;
  
  /* Same again for B */
  for  (  i = 0  ;  i < Nb  &&  B[ i ] <= win[ 1 ]  ;  i++  )  ;
  Nb -= ( Nb  -  i ) ;
  
  
  /*-- Compute STTC --*/
  
  /* All delta-t values */
  for  (  i = 0  ;  i < dtn  ;  i++ )
  {
    
    /* Compute delta-t value */
    dtv = ( double ) i  /  1000.0 ;
    
    /* Save delta-t value */
    if  ( 1  <  nlhs )
      DT[ i ] = dtv ;
    
    /* STTC */
    run_sttc (  &Na  ,  &Nb  ,  &dtv  ,  win  ,  sttc + i  ,  A  ,  B  ) ;
    
  } /* all delta-t */
  
  
  /*-- Return output --*/
  
  /* Return sttc */
  plhs[ STTCARG ] = mxCreateDoubleMatrix (  0  ,  0  ,  mxREAL  ) ;
  mxSetPr (  plhs[ STTCARG ]  ,  sttc  ) ;
  mxSetM  (  plhs[ STTCARG ]  ,   dtn  ) ;
  mxSetN  (  plhs[ STTCARG ]  ,     1  ) ;
  
  /* DT not requested , end now */
  if  ( nlhs  <=  1 )
    return ;
  
  /* Return DT */
  plhs[ DTARG ] = mxCreateDoubleMatrix (  0  ,  0  ,  mxREAL  ) ;
  mxSetPr (  plhs[ DTARG ]  ,   DT  ) ;
  mxSetM  (  plhs[ DTARG ]  ,  dtn  ) ;
  mxSetN  (  plhs[ DTARG ]  ,    1  ) ;
  
  
} /* mexFunction */

