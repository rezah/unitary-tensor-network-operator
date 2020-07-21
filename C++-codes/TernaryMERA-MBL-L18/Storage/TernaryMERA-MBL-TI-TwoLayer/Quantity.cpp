#include"Header.h"
void FullSpectrum(cx_mat * Oascend1,cx_mat *Uni1,cx_mat *Isometry1,cx_mat * Oascend0,cx_mat *Uni0,cx_mat *Isometry0,
ofstream & Energy,double & J,ofstream & SpinOrder,ofstream & SpinOrder1,ofstream &Sectral,ofstream &Sectral1,
vec *Randomvec,ofstream & TimeEvolotion ){
vec Xlvec;
vec Xlvec0;
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
int NumIntChosen;
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
vec s;
cx_mat h;
cx_mat eigvec;
vec eigval;
complex<double> img(0,1);

cx_mat PZ(2,2); PZ.zeros();        PZ(0,0)=1; PZ(1,1)=-1;
cx_mat PI(2,2);  PI.eye(2,2);
cx_mat PX(2,2); PX.zeros();    PX(0,1)=1; PX(1,0)=1;
cx_mat PY(2,2); PY.zeros();    PY(0,1)=-img; PY(1,0)=img;
PX=PX*(1.00/2.00);
PY=PY*(1.00/2.00);
PZ=PZ*(1.00/2.00);



cx_mat * Oout=new cx_mat[1];
cx_mat * Oascend0H=new cx_mat[N0];
cx_mat * Oascend1H=new cx_mat[N0];

cx_mat * Oascend2=new cx_mat[4];


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

NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=0;
 BoundryLeft=0;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];

NumInter=4;
 NumIso=1;
 NumUni=1;
 BoundryRight=0;
 BoundryLeft=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=strans(Oout[0])+Oascend2[0];
cx_mat Oenergy=Oascend2[0];

cx_mat I(Xi1[1],Xi1[1]);  I.eye(Xi1[1],Xi1[1]);
cx_mat I1(Xi1[3],Xi1[3]);  I1.eye(Xi1[3],Xi1[3]);
h=kron(kron(Oascend1[6],kron(I1,I1)),kron(I,I));
h=h+kron(kron(I,Oascend1[1]),kron(I1,kron(I,I)));
h=h+kron(kron(kron(I,I),Oascend1[2]),kron(I,I));
h=h+kron(kron(kron(I,I),I1),kron(Oascend1[3],I));
h=h+kron(kron(kron(I,I),kron(I1,I1)),Oascend1[4]);


h=h+kron(kron(kron(PX,kron(I1,I1)),kron(I1,I1)),PX);
h=h+kron(kron(kron(PY,kron(I1,I1)),kron(I1,I1)),PY);
h=h+kron(kron(kron(PZ,kron(I1,I1)),kron(I1,I1)),PZ);
h=h-Randomvec[0](5)*kron(kron(kron(PI,kron(I1,I1)),kron(I1,I1)),PZ);




PX=PX*(2.00);
PY=PY*(2.00);
PZ=PZ*(2.00);

eig_sym(eigval, eigvec, h,"std");

s.zeros(pow(2,Nreal));
for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2]);i++)
for(int j=0;j<int(Xi1[3]*Xi1[4]*Xi1[2]);j++)
s(i*int(Xi1[3]*Xi1[4]*Xi1[2])+j)=real(Oenergy(i,j));

s=sort(s);

for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2])*int(Xi1[1]*Xi1[2]*Xi1[2]);i++){
Energy<<J;
Energy<<"   "<<s(i)<<"  "<<eigval(i)<<endl;}

cx_mat SigmaZ=kron(kron(kron(kron(PI,PI),kron(PZ,PI)),PI),PI);
cx_mat SigmaZ1=kron(kron(kron(kron(PI,PI),kron(PI,PZ)),PI),PI);

SigmaZ=trans(eigvec)*SigmaZ*eigvec;
SigmaZ1=trans(eigvec)*SigmaZ1*eigvec;

for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2])*int(Xi1[1]*Xi1[2]*Xi1[2]);i++){
SpinOrder<<J;
SpinOrder<<"   "<<eigval(i)<<"  "<<real(SigmaZ(i,i))<<"   "<<real(SigmaZ1(i,i))<<endl;
}

for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2])*int(Xi1[1]*Xi1[2]*Xi1[2]);i++)
for(int j=0;j<int(Xi1[1]*Xi1[2]*Xi1[2])*int(Xi1[1]*Xi1[2]*Xi1[2]);j++){
Sectral<<J;
Sectral<<"   "<<eigval(i)-eigval(j)<<"  "<<norm(real(SigmaZ(i,j)))<<"   "<<norm(real(SigmaZ1(i,j)))<<endl;
}
/*******************************************************************************/
Oascend1H[2]=kron(PZ,PI);
Oascend0H[2]=kron(PI,PZ);

NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=0;
BoundryLeft=0;
NumIntChosen=2;
AscenOne(Oascend1H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend1H[0]=Oout[0];

NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=0;
BoundryLeft=0;
NumIntChosen=2;
AscenOne(Oascend0H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend0H[0]=Oout[0];


for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2]);i++)
for(int j=0;j<int(Xi1[1]*Xi1[2]*Xi1[2]);j++)
{
SpinOrder1<<J;
SpinOrder1<<"   "<<real(Oenergy(i,j))<<"  "<<real(Oascend0H[0](i,j))<<"  "<<real(Oascend1H[0](i,j))<<endl;}

/************************************Spectral***************************************************************************/

Oascend1H[2]=kron(PZ,PI);
Oascend0H[2]=kron(PI,PZ);

NumInter0=1;
NumIso0=0;
NumUni0=0;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend1H,Uni1,Isometry1,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec,Oout,NumIntChosen);
Oascend1H[0]=Oout[0];

NumInter0=1;
NumIso0=0;
NumUni0=0;
BoundryRight0=0;
BoundryLeft0=0;
NumIntChosen=2;
AscenMERAONE(Oascend0H,Uni1,Isometry1,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec,Oout,NumIntChosen);
Oascend0H[0]=Oout[0];




for(int i=0;i<int(Xi1[1]*Xi1[2]*Xi1[2]);i++)
for(int j=0;j<int(Xi1[1]*Xi1[2]*Xi1[2]);j++)
for(int m=0;m<int(Xi1[1]*Xi1[2]*Xi1[2]);m++)
for(int n=0;n<int(Xi1[1]*Xi1[2]*Xi1[2]);n++)
{
Sectral1<<J;
Sectral1<<"   "<<real(Oenergy(i,j))-real(Oenergy(m,n))<<"  "<<
abs(real(Oascend0H[0](m*int(Xi1[1]*Xi1[2]*Xi1[2])+n,i*int(Xi1[1]*Xi1[2]*Xi1[2])+j)))<<
"  "<<abs(real(Oascend1H[0](m*int(Xi1[1]*Xi1[2]*Xi1[2])+n,i*int(Xi1[1]*Xi1[2]*Xi1[2])+j)))<<endl;}
/***************************************************************************************************/
cx_vec Initialvec;
Initialvec.zeros(pow(2,6));
Initialvec=randu<cx_vec>(pow(2,6));
Initialvec=Initialvec/norm(Initialvec);
cout<<norm(Initialvec)<<endl;
cout<<trace(trans(Initialvec)*Initialvec)<<endl;

SigmaZ=kron(kron(kron(kron(PI,PI),kron(PZ,PI)),PI),PI);
SigmaZ1=kron(kron(kron(kron(PI,PI),kron(PI,PZ)),PI),PI);


double Sum;
double Sum2;

for(double i=0.0;i<Time;i=i+TimeSteps)
{
Sum=real(trace(trans(eigvec*diagmat(exp(img*i*eigval))*trans(eigvec)*Initialvec)*SigmaZ*(
eigvec*diagmat(exp(img*i*eigval))*trans(eigvec)*Initialvec)));

Sum2=real(trace(trans(eigvec*diagmat(exp(img*i*eigval))*trans(eigvec)*Initialvec)*SigmaZ1*(
eigvec*diagmat(exp(img*i*eigval))*trans(eigvec)*Initialvec)));

TimeEvolotion<<i<<"  "<<Sum<<"  "<<Sum2<<endl;
}









delete [] Oout;
delete [] Oascend1H;
delete [] Oascend0H;
delete [] Oascend2;
}
