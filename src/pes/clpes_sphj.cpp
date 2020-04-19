/**************************************************************
 *     Purpose: This program computes the spherical Bessel      *
 *              functions jn(x) and jn'(x) using subroutine     *
 *              SPHJ                                            *
 *     Input :  x --- Argument of jn(x)                         *
 *              n --- Order of jn(x)  ( n = 0 to 250 )          * 
 *     Output:  SJ(n) --- jn(x)                                 *
 *              DJ(n) --- jn'(x)                                * 
 *     Example:   x =10.0                                       *
 *                n          jn(x)              jn'(x)          *
 *              --------------------------------------------    *
 *                0    -.5440211109D-01    -.7846694180D-01     * 
 *                1     .7846694180D-01    -.7009549945D-01     *
 *                2     .7794219363D-01     .5508428371D-01     *
 *                3    -.3949584498D-01     .9374053162D-01     *
 *                4    -.1055892851D+00     .1329879757D-01     *
 *                5    -.5553451162D-01    -.7226857814D-01     *
 * ------------------------------------------------------------ *
 * REFERENCE: "Fortran Routines for Computation of Special      *
 *             Functions,                                       *
 *             jin.ece.uiuc.edu/routines/routines.html".        *
 *                                                              *
 *                     Visual C++ Release By J-P Moreau, Paris. *
 *                                 (www.jpmoreau.fr)            *
 ***************************************************************/
/*
*/


void clpes::SPHJ(int N, double X, int *NM, double *SJ) {
  /*     =======================================================
!      Purpose: Compute spherical Bessel functions jn(x) and
!               their derivatives
!      Input :  x --- Argument of jn(x)
!               n --- Order of jn(x)  ( n = 0,1,... )
!      Output:  SJ(n) --- jn(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      ======================================================= */
  int K, M; double CS, F, F0, F1, SA, SB;
        
  *NM=N;
  if (fabs(X) < 1e-100) {
    for (K=0; K<=N; K++) {
      SJ[K]=0.0;
    } 
    SJ[0]=1.0;
    return;
  }
  SJ[0]=sin(X)/X;
  SJ[1]=(SJ[0]-cos(X))/X;
  if (N >= 2) {
    SA=SJ[0];
    SB=SJ[1];
    M=MSTA1(X,200);
    if (M < N)
      *NM=M;
    else
      M=MSTA2(X,N,15);
    F0=0.0;
    F1=1.0-100;
    for (K=M; K>-1; K--) {
      F=(2.0*K+3.0)*F1/X-F0;
      if (K <= *NM)  SJ[K]=F;
      F0=F1;
      F1=F;
    }
    if (fabs(SA) > fabs(SB))  CS=SA/F;
    if (fabs(SA) <= fabs(SB)) CS=SB/F0;
    for (K=0; K<=*NM; K++)  SJ[K] *= CS;
  }      
}

int clpes::MSTA1(double X, int MP) {

  /*     ===================================================
!      Purpose: Determine the starting point for backward  
!               recurrence such that the magnitude of    
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point   
!      =================================================== */

  double A0,F,F0,F1; int IT,NN,N0,N1;
  A0=fabs(X);
  N0=floor(1.1*A0)+1;
  F0=ENVJ(N0,A0)-MP;
  N1=N0+5;
  F1=ENVJ(N1,A0)-MP;
  for (IT=1; IT<=20; IT++) {             
    NN=N1-(N1-N0)/(1.0-F0/F1);                  
    F=ENVJ(NN,A0)-MP;
    if (abs(NN-N1) < 1) goto e20;
    N0=N1;
    F0=F1;
    N1=NN;
    F1=F;
  }
 e20:return NN;
}

int clpes::MSTA2(double X, int N, int MP) {

  /*     ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      =================================================== */

  double A0,EJN,F,F0,F1,HMP,OBJ; int IT,N0,N1,NN;
  A0=fabs(X);
  HMP=0.5*MP;
  EJN=ENVJ(N,A0);
  if (EJN <= HMP) {
    OBJ=MP;
    N0=floor(1.1*A0);
  }
  else {
    OBJ=HMP+EJN;
    N0=N;
  }
  F0=ENVJ(N0,A0)-OBJ;
  N1=N0+5;
  F1=ENVJ(N1,A0)-OBJ;
  for (IT=1; IT<=20; IT++) {
    NN=N1-(N1-N0)/(1.0-F0/F1);
    F=ENVJ(NN,A0)-OBJ;
    if (abs(NN-N1) < 1) goto e20;
    N0=N1;
    F0=F1;
    N1=NN;
    F1=F;
  }
 e20:    return NN+10;
}


double clpes::ENVJ(int N, double X) {
  return (0.5*log(6.28*N)-N*log(1.36*X/N));
}
