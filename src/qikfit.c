/* qikfit subroutine  
from Andrew Martin's BioPlib package
http://www.bioinf.org.uk/software/profit/index.html
below is original header:

*************************************************************************

   Program:    
   File:       fit.c
   
   Version:    V1.4R
   Date:       03.06.97
   Function:   Perform least squares fitting of coordinate sets
   
   Copyright:  (c) SciTech Software 1993-7
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead, Surrey, KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      amartin@stagleys.demon.co.uk
               martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   This code performs least squares fitting of two coordinate set using
   the method of McLachlan as modified by Sutcliffe.

**************************************************************************
 to adapt code to pdb2pt:
  * variable type for "column" changed from bool to int
     (pass '1' instead of TRUE for pdb2icos matrix element order)
  * variable type REAL changed to double, eliminating mathType.h include
  * macro ABS replaced with math function fabs, eliminating macros.h include
  C Lawson August 2006
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10


void qikfit (double umat[3][3], double rm[3][3], int column) {
   
   double  rot[3][3], turmat[3][3],
         c[3][3], coup[3], dir[3], step[3], v[3],
         rtsum,rtsump, rsum,
         stp,stcoup,
         ud,tr,ta,cs,sn,ac,
         delta,deltap, gfac,
         cle,clep;
   int   i,j,k,l,m,
         jmax, ncyc, nsteep, nrem;

   /* Rotate repeatedly to reduce couple about initial direction to zero.
      Clear the rotation matrix
   */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp)                             */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step
      */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple                                                         */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle     = sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup                                   */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step                              */
      rtsump = rtsum;
      deltap = delta;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous                             */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step                                                */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction                           */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction                              */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then position is unstable, so don't stop               */
      if((cs>0.0) && (fabs(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part
      */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part                                                */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix                                   */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor                                             */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop                      */
      if(fabs(delta)<SMALL) break;

      /* Next cycle                                                     */
   }

   rsum = rtsum;

   /* Copy rotation matrix for output                                   */
   if(column==1)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}
