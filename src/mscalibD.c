#include "mscalibD.h"


void getaccD(
  double *pl
  ,Sint *lpll
  ,double *cal
  ,Sint *lcall
  ,double *error
  ,int *plind
  ,int *calind
  ,Sint *ind
  ,Sint *ppm
  )
{
  	/*pl peaklist 1*/
  	/*lpl length peaklist1*/
  	/*cal peaklist 2*/
  	/*lcal length of peaklist 2*/
  	/*mods modifications*/
  	/*lmods length modifications*/
  	/*error - allowed error*/
  	/*plind - return array with peaklist indices*/
  	/*calind -return array with calib indices*/
  	/*mind - return array with modification idices*/
  	/*length of return arrays*/
  	/*ppm if 1 than use ppm else normal error*/
    int lpl=*lpll;
    int lcal=*lcall;
    int nind = *ind;
    double terror = *error;
  	int p = 0, c=0;
  	for(p=0;p<lpl;p++)
  	{
  	  for(c=0;c<lcal;c++)
  	  {
  	    if(*ppm)
  	    {
  	      if((pl[p] + (terror*1000000)/pl[p]) < cal[c]){break;} /*cal[c] is to large therfore no more matches*/
          while((pl[p]-(terror*1000000)/pl[p]) > cal[c] && c < (lcal-1)){c++;}/*cal[c] is to small therefore no more matches */     
  	      if( ( fabs(pl[p]- cal[c])*1000000)/pl[p] < terror)
	        {
	            plind[nind]= p;
	            calind[nind]= c;
	            nind = nind+1;
             /* printf("HIT: %i %i %i :  %f %f %f %f\n",p,c,nind, pl[p], cal[c],fabs(pl[p]-cal[c]), terror);*/
	        } 
	       }
  	    else
  	    {
  	      if((pl[p] + terror) < cal[c]){break;} /*increases pl[p]*/
          while((pl[p]-terror) > cal[c] && c < (lcal-1) ){c++;}
          if( fabs(pl[p]-cal[c])<terror)
  	      {
  	          plind[nind]  = p;
	            calind[nind] = c;
	            nind  =   nind + 1;
	           /* printf("HIT: %i %i %i :  %f %f %f %f\n",p,c,nind, pl[p], cal[c],fabs(pl[p] - cal[c]), terror);*/
  	      }
  	    }
  	  }
  	}
  	*ind=nind;
}



void getaccU(
  double *pl
  ,Sint *lpll
  ,double *cal
  ,Sint *lcall
  ,double *error
  ,int *plind
  ,int *calind
  ,Sint *ind
  ,Sint *ppm
  )
{
        int p=0,c=0;
        int hits=*ind;
        double window=*error;
        int lpl = *lpll;
        int lcal= *lcall;
        
        
        if ( window < 0 ) { printf("FATAL ERROR: Need a positive shift value.\n"); exit(1); }

  if(*ppm)
  {
    while ( p<lpl && c<lcal ) {
        if ( pl[p] >= ( cal[c]-(window*cal[c])/1000000 ) && pl[p] <= ( cal[c]+(window*cal[c])/1000000 ) ) { /* did we hit something ??? */
                plind[hits]=p;
                calind[hits]=c;
                /*printf("HIT: %i %i %i :  %f %f %f %f\n",p,c,hits, pl[p], cal[c],pl[p] - cal[c], 0.6);*/
                hits++; p++; c++; 
        } else { 
          while (   pl[p] <  ( cal[c]-(window*cal[c])/1000000 ) && p < lpl ) { p++; }   /* a<b -> move a */
          while ( ( pl[p]-(window*cal[c])/1000000) >  cal[c]   && c < lcal ) { c++; }
        }
    }
    
  }else{      
     while ( p<lpl && c<lcal ) {
     
      if ( pl[p] >= ( cal[c] - window ) && pl[p] <= ( cal[c]+window ) ) { /* did we hit something ???*/
              plind[hits]=p;
              calind[hits]=c;
              /*printf("HIT: %i %i %i :  %f %f %f %f\n", p , c , hits ,  pl[p] , cal[c], pl[p] - cal[c] , 0.6);*/
              hits++; p++; c++; 
      } else { 
        while (   pl[p] <  ( cal[c]-window ) && p < lpl ) { p++; }   /* a<b -> move a*/
        while ( ( pl[p]-window ) >  cal[c]   && c < lcal ) { c++; }
      }
    }
  }
  *ind=hits;
} /* match_peaks*/
