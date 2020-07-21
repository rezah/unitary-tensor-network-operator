#include"Header.h"
using namespace std;
int main()
{
cx_mat * Oascend0=new cx_mat[7];
cx_mat * Oascend1=new cx_mat[7];
cx_mat * Oascend1H=new cx_mat[7];
cx_mat * Oascend2=new cx_mat[2];
cx_mat * Isometry0=new cx_mat[int(N0/3)];
cx_mat * Isometry1=new cx_mat[3];
cx_mat * Uni0=new cx_mat[int(N0/3)-1];
cx_mat * Uni1=new cx_mat[2];
cx_mat * Rho0=new cx_mat[N0];
cx_mat * Rho1=new cx_mat[int(N0/3)-1];
cx_mat * Rho2=new cx_mat[1];
cx_mat * Env=new cx_mat[1];
cx_mat * EnvU=new cx_mat[4];
cx_mat * EnvIso=new cx_mat[2];
cx_mat * Oout=new cx_mat[1];
vec * CouplingRandoms=new vec[1];
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


cx_mat U;
cx_mat V;
vec s;
cx_mat h;
cx_mat eigvec;
vec eigval;
double Gamma=1;
double E1,E2;
double NormReal;
cx_mat DX;
int Count;
int BreakWhile=1;
cx_mat Unih;
cx_mat IsoLh;
cx_mat IsoRh;
cx_mat Z;
double J=CouplingJ;
double traceH2;

ofstream Energy;
Energy.open("Files/Energy.text");

ofstream File;
File.open("Files/File.text");

ofstream AccuracyFile;
AccuracyFile.open("Files/AccuracyFile.text");

ofstream SpinOrder;
SpinOrder.open("Files/SpinOrder.text");

ofstream SpinOrder1;
SpinOrder1.open("Files/SpinOrder1.text");


ofstream Sectral;
Sectral.open("Files/Sectral.text");


ofstream Sectral1;
Sectral1.open("Files/Sectral1.text");


ofstream ConfEnergy;
ConfEnergy.open("Files/ConfEnergy.text");

ofstream TimeEvolotion;
TimeEvolotion.open("Files/TimeEvolotion.text");




vec * Checker=new vec[1];
Checker[0].zeros(4);
Checker[0](1)=2200.594;

for(int x0=0;x0<Points;x0++){


if(x0==0)
;
else
J=J+(CSIGN*Grid);

for(int x1=0;x1<Realization;x1++){

UniInitial(Uni0,Uni1,Isometry0,Isometry1);
OInitial(Oascend0,Oascend1,Uni0,Isometry0,J,traceH2,ConfEnergy,CouplingRandoms);

w=1;
for(int q=0;q<=w;q++){

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
cout<<setprecision(16)<<"error="<<traceH2-abs((trace(Oascend2[0]*strans(Oascend2[0])))*(1.00/pow(2,Nreal)))<<endl;

cout<<setprecision(16)<<"ET="<<abs((trace(Oascend2[0]*strans(Oascend2[0])))*(1.00/pow(2,Nreal)))<<endl;


AccuracyFile<<q<<"  "<<traceH2-abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;

OptimizeUp(Oascend0,Uni0,Isometry0,Oascend1,Uni1,Isometry1,Checker,Oascend2);


if(abs(Checker[0](0))>0.000001){
if( (((abs(Checker[0](0)-Checker[0](1)))/(abs(Checker[0](0))))< Accuracy) && (q==w) )
;
else if(q==w)
w=w+1;
}else{
if( (abs(Checker[0](0)-Checker[0](1))< 0.00000000001) && (q==w) )
;
else if(q==w)
w=w+1;
}

}


cout<<endl<<"TotalIteration="<<w<<endl;
File<<w<<endl;
FullSpectrum(Oascend1,Uni1,Isometry1,Oascend0,Uni0,Isometry0,Energy,J,SpinOrder,SpinOrder1,Sectral,Sectral1,CouplingRandoms,
TimeEvolotion);

}
}


delete [] Oascend0;
delete [] Oascend1;
delete [] Oascend1H;
delete [] Oascend2;
delete [] Isometry0;
delete [] Isometry1;
delete [] Uni0;
delete [] Uni1;
delete [] Rho0;
delete [] Rho1;
delete [] Rho2;
delete [] Env;
delete [] EnvU;
delete [] EnvIso;
delete [] Oout;
delete [] Checker;
delete [] CouplingRandoms;
return 0;
}





