void vypocet();
int compare_EPS(double X[], int n);
void partial();

const double EPS=1e-20;       // presnost vypoctu
const double h=0.1;           // krok reseni
const double T0=0.0;          // pocatecni cas
const double TMAX=1.0;        // max cas
#define MAX_ORD 20            // maximalni rad metody pri vypoctu
#define N 11                  // pocet reseni TMAX / h
double y[N]={0,},z[N]={0,},u[N]={0,}, v[N]={0,};        // reseni
double DY[MAX_ORD]={0,}, DZ[MAX_ORD]={0,}, DU[MAX_ORD]={0},    // prirustky
       DV[MAX_ORD]={0,};
const unsigned short int num_less=3;  // pocet hodnot mensich EPS pro ukonceni vyp.

double vysl_DY=0,vysl_DZ=0,vysl_DU=0,vysl_DV=0;


// konec div.h