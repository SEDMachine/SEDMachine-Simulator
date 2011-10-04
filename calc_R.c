#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"

#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

main()
{
  float *x,*y,R;
  FILE *fin;
  int i,istart,N=55;
  int status = 0, ii, jj;
  long  fpixel = 1, naxis = 2, nelements, exposure;
  long naxes[2] = { 500, 20 };   /* image is 300 pixels wide by 200 rows */
  short array[20][500];
  fitsfile *fptr;
  
  x=(float *)calloc(N,sizeof(float));
  y=(float *)calloc(N,sizeof(float));
  

  fin=fopen("Rtest.dat","r");
  for(i=0;i<N;i++)
    fscanf(fin,"%f %f",&x[i],&y[i]);
  fclose(fin);

  for(i=0;i<(N);i++) {
    R = x[i]/(x[i]-x[i-1]);
    printf("%2.3f\t%2.3f\t%2.3f\t%d\n",x[i],y[i]*1000.0+1500,y[i]*1000.0/13.5+120,round(y[i]*1000.0/13.5+120));
  }
  

  status = 0;         /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, "testfile.fits", &status);   /* create new file */
  
  /* Create the primary array image (16-bit short integer pixels */
  fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
  
  /* Write a keyword; must pass the ADDRESS of the value */
  //exposure = 1500.;
  //fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
  //	  "Total Exposure Time", &status);
  
  /* Initialize the values in the image with a linear ramp function */
  for (jj = 0; jj < naxes[1]; jj++)
    for (ii = 0; ii < naxes[0]; ii++)
      array[jj][ii] = 0;
  

  /* put the values in */
  for (jj = 9; jj < 13; jj++) {
    //istart=0;
    for (ii = 0; ii < naxes[0]; ii++) {
      for(i=0;i<(N-1);i++) {
	if(ii>=round(y[i]*1000.0/13.5+120) && ii<round(y[i+1]*1000.0/13.5+120)) {
	  array[jj][ii] = x[i]*1000.0;
	  //istart++;
	}
	//else
	//break;

      }
    }
  }
  nelements = naxes[0] * naxes[1];          /* number of pixels to write */
  
  /* Write the array of integers to the image */
  fits_write_img(fptr, TSHORT, fpixel, nelements, array[0], &status);
  
  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */
  return( status );
  
  free(x); free(y);
}
  
