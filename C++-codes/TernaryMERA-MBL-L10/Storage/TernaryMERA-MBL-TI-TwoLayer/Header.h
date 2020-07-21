#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED
#include <iostream>
#include<complex>
#include<math.h>
#include <armadillo>
#include<fstream>
#include <iostream>
#include <iomanip>
#include <random>
using namespace std;
using namespace arma;
const int N0=6;
const int Nreal=6;

const int Xi0[6]={2,2,2,2,2,2}; /* Last two elements must be equal*/
const int Xi1[6]={2,2,2,2,2,2}; /* {x0,x1,x2,x2b,x3,x4}*/
const int Xi2[6]={1,1,4,2,2,2}; /* {x0,x1,x2,x2b,x3,x4}*/
const int Xi3[6]={2,2,2,4,1,1}; /* {x0,x1,x2,x2b,x3,x4}*/

const int Xi4[6]={1,1,1,1,1,4}; /* {x0,x1,x2,x2b,x3,x4}*/
const int Xi5[6]={4,1,1,1,1,1}; /* {x0,x1,x2,x2b,x3,x4}*/

const int Points=1;
const double Grid=0.0;
const double CouplingJ=2.0;
const double CSIGN=1;
const int Realization=1;

const int TotallRepeat=500;
const int NumIterationUni0=5;
const int NumIterationIso0=5;
const int NumIterationUP=10;


const int UseStorage=0;
const int NewGeneration=1;
const double Accuracy=1.0e-9;
const double TimeSteps=1.0e-2;
const double Time=1.0e3;





void FullSpectrum(cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *, ofstream &,double &,ofstream &,ofstream &,
ofstream &,ofstream &,vec *,ofstream & );
void OptimizeU(cx_mat*,cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *);
void OptimizeIso(cx_mat*,cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *);
void OptimizeUp(cx_mat*,cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,vec*,cx_mat*);
void OInitial(cx_mat*,cx_mat*,cx_mat*,cx_mat*,double&,double &,ofstream &,vec *);
void UniInitial(cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void Ascen(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,vec,cx_mat *);
void AscenOne(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,int,vec,cx_mat *);
void AscendTotal(cx_mat *,cx_mat*,cx_mat*,cx_mat*);
void EnvUni(cx_mat*,cx_mat*,cx_mat*,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoLeft(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoLeftMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoRight(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);
void AscenMERA(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,vec,cx_mat *);
void AscenMERAONE(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,vec,cx_mat *,int);
void Dscen(cx_mat *,cx_mat*,cx_mat*,cx_mat* ,int,int,int,int,int,int,int,int,vec);
void Dscenone(cx_mat *,cx_mat*,cx_mat*,cx_mat* ,int,int,int,int,int,int,int,vec);
void RhoLanczos(cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void EnvIsoRightMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,vec,cx_mat *);
void EnvUniMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,vec,cx_mat *);
void OptUni(cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,int);
void OptIso(cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,int);
void OptTop(cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void OptTopFinite(cx_mat*,cx_mat*,cx_mat*,cx_mat*,double&);
#endif // HEADER_H_INCLUDED
