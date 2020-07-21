#include"Header.h"
void FullSpectrum(cx_mat * Oascend1,cx_mat *Uni1,cx_mat *Isometry1,cx_mat * Oascend0,cx_mat *Uni0,cx_mat *Isometry0,
ofstream & Energy,double & J,ofstream & SpinOrder,ofstream & SpinOrder1,ofstream &Sectral,ofstream &Sectral1 ){

int NumIntChosen;

int NumInter;
int NumIso;
int NumUni;
int BoundryRight;
int BoundryLeft;
int NumInterLeft;
int NumInterRight;
int NumRhoup;
int NumRhodown;
int NumRhoLeft;
int NumRhoRight;
int w;

int NumInter0;
int NumIso0;
int NumUni0;
int BoundryRight0;
int BoundryLeft0;
int NumInterLeft0;
int NumInterRight0;
int NumRhoup0;
int NumRhodown0;
int NumRhoLeft0;
int NumRhoRight0;

complex<double> img(0,1);

cx_mat PZ(2,2); PZ.zeros();        PZ(0,0)=1; PZ(1,1)=-1;
cx_mat PI(2,2);  PI.eye(2,2);
cx_mat PX(2,2); PX.zeros();    PX(0,1)=1; PX(1,0)=1;
cx_mat PY(2,2); PY.zeros();    PY(0,1)=-img; PY(1,0)=img;


vec Xlvec;
vec Xlvec0;

Xlvec0.zeros(6);
Xlvec0(0)=Xi0[0];
Xlvec0(1)=Xi0[1];
Xlvec0(2)=Xi0[2];
Xlvec0(3)=Xi0[3];
Xlvec0(4)=Xi0[4];
Xlvec0(5)=Xi0[5];

Xlvec.zeros(6);
Xlvec(0)=Xi1[0];
Xlvec(1)=Xi1[1];
Xlvec(2)=Xi1[2];
Xlvec(3)=Xi1[3];
Xlvec(4)=Xi1[4];
Xlvec(5)=Xi1[5];

vec s;
cx_mat * Oascend0H=new cx_mat[N0];
cx_mat * Oascend1H=new cx_mat[N0];
cx_mat * Oout=new cx_mat[1];


AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=1;
BoundryLeft=1;
NumInterLeft=0;
NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cx_mat Oenergy=Oout[0];



s.zeros(pow(2,N0));
for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0]);i++)
for(int j=0;j<int(Xi1[3]*Xi1[4]*Xi1[5]);j++)
s(i*int(Xi1[3]*Xi1[4]*Xi1[5])+j)=real(Oout[0](i,j));
s=sort(s);
for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0])*int(Xi1[3]*Xi1[4]*Xi1[5]);i++){
Energy<<J;
Energy<<"   "<<s(i)<<endl;}
/*********************************************/

Oascend0H[8]=kron(PZ,PI);
NumInter0=7;
NumIso0=2;
NumUni0=2;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend0H,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,
NumInterRight0,Xlvec0,Oout,NumIntChosen);
Oascend1H[2]=Oout[0];


Oascend0H[8]=kron(PI,PZ);
NumInter0=7;
NumIso0=2;
NumUni0=2;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend0H,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,
NumInterRight0,Xlvec0,Oout,NumIntChosen);
Oascend0H[2]=Oout[0];


NumInter=1;
NumIso=0;
NumUni=0;
NumIntChosen=2;
AscenOne(Oascend1H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
NumIntChosen,Xlvec,Oout);
Oascend1H[0]=Oout[0];


NumInter=1;
NumIso=0;
NumUni=0;
NumIntChosen=2;
AscenOne(Oascend0H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
NumIntChosen,Xlvec,Oout);
Oascend0H[0]=Oout[0];




for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0]);i++)
for(int j=0;j<int(Xi1[3]*Xi1[4]*Xi1[5]);j++)
{
SpinOrder<<J;
SpinOrder<<"   "<<real(Oenergy(i,j))<<"  "<<real(Oascend1H[0](i,j))<<"  "<<real(Oascend0H[0](i,j))<<endl;}

/****************************************/

UniInitial(Uni0,Uni1,Isometry0,Isometry1);


AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=1;
BoundryLeft=1;
NumInterLeft=0;
NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oenergy=Oout[0];



s.zeros(pow(2,N0));
for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0]);i++)
for(int j=0;j<int(Xi1[3]*Xi1[4]*Xi1[5]);j++)
s(i*int(Xi1[3]*Xi1[4]*Xi1[5])+j)=real(Oout[0](i,j));
s=sort(s);
for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0])*int(Xi1[3]*Xi1[4]*Xi1[5]);i++){
Sectral<<J;
Sectral<<"   "<<s(i)<<endl;}


Oascend0H[8]=kron(PZ,PI);
NumInter0=7;
NumIso0=2;
NumUni0=2;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend0H,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,
NumInterRight0,Xlvec0,Oout,NumIntChosen);
Oascend1H[2]=Oout[0];


Oascend0H[8]=kron(PI,PZ);
NumInter0=7;
NumIso0=2;
NumUni0=2;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend0H,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,
NumInterRight0,Xlvec0,Oout,NumIntChosen);
Oascend0H[2]=Oout[0];


NumInter=1;
NumIso=0;
NumUni=0;
NumIntChosen=2;
AscenOne(Oascend1H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
NumIntChosen,Xlvec,Oout);
Oascend1H[0]=Oout[0];


NumInter=1;
NumIso=0;
NumUni=0;
NumIntChosen=2;
AscenOne(Oascend0H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
NumIntChosen,Xlvec,Oout);
Oascend0H[0]=Oout[0];




for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[0]);i++)
for(int j=0;j<int(Xi1[3]*Xi1[4]*Xi1[5]);j++)
{
SpinOrder1<<J;
SpinOrder1<<"   "<<real(Oenergy(i,j))<<"  "<<real(Oascend1H[0](i,j))<<"  "<<real(Oascend0H[0](i,j))<<endl;}


delete [] Oout;
delete [] Oascend1H;
delete [] Oascend0H;
}

