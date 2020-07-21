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
const int N0=10;
const int Xi0[6]={1,1,2,2,1,1}; /* Last two elements must be equal*/
const int Xi1[6]={2,2,8,8,2,2}; /* {x0,x1,x2,x2b,x3,x4}*/

const int Points=1;
const double Grid=0.10;
const double CouplingJ=5.0;
const double CSIGN=1;   //if it's equal to zero, fidelity is not being calculated.
const int TotallRepeat=0;

const int UseStorage=0;
const int NewGeneration=1;
const double Accuracy=1.0e-4;
const double AccuracyTop=1.0e-8;

const double DeltaFidelity=0.000;   //if it's equal to zero, fidelity is not being calculated.
const int PointsFidelity=1;   //>=2 otherwise must be always equal to 1
const double gridFidelity=0.0000;

const double JY=0;
const double JZ=0;
void GSF( cx_mat * , cx_mat * , ofstream &,ofstream&,ofstream&,double&);
void GSFMORE( cx_mat * , cx_mat * , ofstream &,ofstream&,ofstream&,double&);
void Storing(cx_mat* ,cx_mat* );
void Quantity(ofstream&,ofstream&,ofstream&,ofstream&,ofstream&,ofstream&
,ofstream&,cx_mat*,cx_mat*,cx_mat*,double ,int,double&,ofstream& );
void OInitial(cx_mat*,cx_mat*,double&,double & );
void UniInitial(cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void Ascen(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,vec,cx_mat *);
void EnvUni(cx_mat*,cx_mat*,cx_mat*,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoLeft(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoRight(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,int,int,vec,cx_mat *);

void AscenMERA(cx_mat*,cx_mat*,cx_mat*,int ,int,int,int,int,int,int,vec,cx_mat *);

void Dscen(cx_mat *,cx_mat*,cx_mat*,cx_mat* ,int,int,int,int,int,int,int,int,vec);
void RhoLanczos(cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void EnvIsoLeftMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,vec,cx_mat *);
void EnvIsoRightMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,int,int,vec,cx_mat *);
void EnvUniMERA(cx_mat*,cx_mat*,cx_mat* ,cx_mat*,int,int,int,int,vec,cx_mat *);
void OptUni(cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,int);
void OptIso(cx_mat *,cx_mat *,cx_mat *,cx_mat *,cx_mat *,int);
void OptTop(cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*);
void OptTopFinite(cx_mat*,cx_mat*,cx_mat*,cx_mat*,double&);
#endif // HEADER_H_INCLUDED
