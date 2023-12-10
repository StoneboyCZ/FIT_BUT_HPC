#include <stdio.h>
#include <math.h>

#include "div.h"
//*********************************************
// je treba nejprve nastavit y[0],z[0],u[0];
//*********************************************

//----------------------------------------------------------------------
// kontrola ukonceni vypoctu
// num_less hodnot mensich nez EPS
// 1 pokud tri posledni hodnoty mensi EPS
int compare_EPS(double X[], int n)
{
//  printf("Test EPS :n=%d\n",n);
  int res=1,i=0;
  if (n<num_less) return 0;
  while (i<num_less && res)
  {
     res = fabs(X[n-i])<EPS;
     i++;
  }
  return res;
}

//-------------------------------------------------------------------
// hlavni vypocet
void vypocet()
{
  int i=0;    // ridici promenna pro vysledky y,z,u
  double t=T0;
  while (t<TMAX-h)
  {
//     printf("DY[0]=y[%d]=%g\n",i,y[i]);
     DY[0]=y[i];
     DZ[0]=z[i];
     DU[0]=u[i];
     DV[0]=v[i];
     partial();
     y[i+1] = y[i] + vysl_DY;
     z[i+1] = z[i] + vysl_DZ;
     u[i+1] = u[i] + vysl_DU;
     v[i+1] = v[i] + vysl_DV;
     i++;
     for (int pom=0;pom<MAX_ORD;pom++)
     {
       DY[pom]=0;
       DZ[pom]=0;
       DU[pom]=0;
       DV[pom]=0;
     }
     vysl_DY=0;
     vysl_DZ=0;
     vysl_DU=0;
     vysl_DV=0;
     t+=h;
  }
}

//------------------------------------------------------------------
// rekurzni generovani koeficientu DY = h * DU * DZ
void generate(double *result,int i,int j)
{
  if (i>=1)
  {
    *result+=DU[i]*DY[j];
    printf("result: %g, i: %d, j: %d, du: %g, dy: %g\n", *result, i, j, DU[i], DY[j]);
    i--;
    j++;
    generate(result,i,j);
  }
}

//------------------------------------------------------------------
// vypocet DY,DZ, DU
void partial()
{
  int i=0;
  double result;
  while (i<MAX_ORD && !compare_EPS(DY,i) && !compare_EPS(DZ,i)
	 && !compare_EPS(DU,i))
  {
      result=0.0;
      DU[i+1]=(h/(double)(i+1))*(DU[i]);  // e^t
      DV[i+1]=(h/(double)(i+1))*(DZ[i]);   // -sin t
      DZ[i+1]=(h/(double)(i+1))*(-DV[i]);  // cos t
      generate(&result,i,0);
      printf("%d\n",i);
      DY[i+1]=(h/(double)((i+1)*DU[0]))*(DZ[i]-result);
      i++;
      vysl_DY+=DY[i];
      vysl_DU+=DU[i];       // e^-i
      vysl_DZ+=DZ[i];       // cos t
      vysl_DV+=DV[i];
  }
}

//----------------------------------------------------------------
//----------------------------------------------------------------
//------------------- M A I N ------------------------------------

int main()
{
 // zadani: y'= cos t / e^t y(0)=1; dt=0.1, TMAX=1
 //         y'=  z  /  u
 
 y[0]=1;
 u[0]=1;
 z[0]=1;
 v[0]=0;   // nastaveni pocatecnich podminek

 vypocet();
 
 printf("t                y                u              z\n");

 for (int i=0; i<N;i++)
   printf("%2d         %.10g     %.10g    %.10g\n",i,y[i],u[i],z[i]);

 printf("=========================================================\n");

 return 0;
}

// konec div.cpp