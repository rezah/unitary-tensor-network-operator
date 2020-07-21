#include"Header.h"
using namespace std;
int main()
{
cx_mat * Oascend0=new cx_mat[N0];
cx_mat * Oascend1=new cx_mat[5];
cx_mat * Oascend1H=new cx_mat[5];
cx_mat * Oascend2=new cx_mat[2];

cx_mat * Isometry0=new cx_mat[int(N0/3)];
cx_mat * Isometry1=new cx_mat[2];

cx_mat * Uni0=new cx_mat[int(N0/3)-1];
cx_mat * Uni1=new cx_mat[1];

cx_mat * Rho0=new cx_mat[N0];
cx_mat * Rho1=new cx_mat[int(N0/3)-1];
cx_mat * Rho2=new cx_mat[1];


cx_mat * Env=new cx_mat[1];
cx_mat * EnvU=new cx_mat[4];
cx_mat * EnvIso=new cx_mat[2];
cx_mat * Oout=new cx_mat[1];

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
vec RealizationError;

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

ofstream Error;
Error.open("Files/Error.text");


vec * Checker=new vec[1];
Checker[0].zeros(4);
Checker[0](1)=2200.594;

for(int x0=0;x0<Points;x0++){


if(x0==0)
;
else
J=J+(CSIGN*Grid);

RealizationError.zeros(Realization);
for(int x1=0;x1<Realization;x1++){

UniInitial(Uni0,Uni1,Isometry0,Isometry1);
OInitial(Oascend0,Oascend1,J,traceH2,ConfEnergy);
if(UseStorage==1)
{
Isometry0[0].load("Storage/Isometry0[0].mat",arma_binary);
Isometry0[1].load("Storage/Isometry0[1].mat",arma_binary);
Isometry0[2].load("Storage/Isometry0[2].mat",arma_binary);
Isometry0[3].load("Storage/Isometry0[3].mat",arma_binary);
Isometry0[4].load("Storage/Isometry0[4].mat",arma_binary);
Isometry0[5].load("Storage/Isometry0[5].mat",arma_binary);
Isometry1[0].load("Storage/Isometry1[0].mat",arma_binary);
Isometry1[1].load("Storage/Isometry1[1].mat",arma_binary);
Uni0[0].load("Storage/Uni0[0].mat",arma_binary);
Uni0[1].load("Storage/Uni0[1].mat",arma_binary);
Uni0[2].load("Storage/Uni0[2].mat",arma_binary);
Uni0[3].load("Storage/Uni0[3].mat",arma_binary);
Uni0[4].load("Storage/Uni0[4].mat",arma_binary);
Uni1[0].load("Storage/Uni1[0].mat",arma_binary);
}


AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);

complex<double> Htopowerto=trace(Oascend1[0]*Oascend1[0])*pow(8,4);
Htopowerto=Htopowerto+trace(Oascend1[1]*Oascend1[1])*pow(8,4);
Htopowerto=Htopowerto+trace(Oascend1[2]*Oascend1[2])*pow(8,4);
Htopowerto=Htopowerto+trace(Oascend1[3]*Oascend1[3])*pow(8,4);
Htopowerto=Htopowerto+trace(Oascend1[4]*Oascend1[4])*pow(8,4);
cout<<endl<<real(Htopowerto*(1.00/pow(2,N0)))<<endl;
traceH2=real(Htopowerto*(1.00/pow(2,N0)));

w=1;
for(int q=0;q<=w;q++){


AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<setprecision(16)<<"ET="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0)))<<endl;


cout<<"Uni"<<endl;
OptimizeU(Oascend0,Uni0,Isometry0,Oascend1,Uni1,Isometry1);
AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<setprecision(16)<<"ET="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0)))<<endl;
cout<<"Iso"<<endl;
OptimizeIso(Oascend0,Uni0,Isometry0,Oascend1,Uni1,Isometry1);
AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<setprecision(16)<<"ET="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0)))<<endl;
cout<<setprecision(16)<<"error="<<traceH2-abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0)))<<endl;
RealizationError(x1)=abs(traceH2-abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0))));
AccuracyFile<<q<<"  "<<abs(traceH2-abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0))))<<endl;




cout<<"Up"<<endl;
OptimizeUp(Oascend0,Uni0,Isometry0,Oascend1,Uni1,Isometry1,Checker);

AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<endl<<setprecision(16)<<"error="<<traceH2-abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,N0)))<<endl;

if((q==w)  || (q%30==0)){
Isometry0[0].save("Storage/Isometry0[0].mat",arma_binary);
Isometry0[1].save("Storage/Isometry0[1].mat",arma_binary);
Isometry0[2].save("Storage/Isometry0[2].mat",arma_binary);
Isometry0[3].save("Storage/Isometry0[3].mat",arma_binary);
Isometry0[4].save("Storage/Isometry0[4].mat",arma_binary);
Isometry0[5].save("Storage/Isometry0[5].mat",arma_binary);
Isometry1[0].save("Storage/Isometry1[0].mat",arma_binary);
Isometry1[1].save("Storage/Isometry1[1].mat",arma_binary);
Uni0[0].save("Storage/Uni0[0].mat",arma_binary);
Uni0[1].save("Storage/Uni0[1].mat",arma_binary);
Uni0[2].save("Storage/Uni0[2].mat",arma_binary);
Uni0[3].save("Storage/Uni0[3].mat",arma_binary);
Uni0[4].save("Storage/Uni0[4].mat",arma_binary);
Uni1[0].save("Storage/Uni1[0].mat",arma_binary);
}

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
FullSpectrum(Oascend1,Uni1,Isometry1,Oascend0,Uni0,Isometry0,Energy,J,SpinOrder,SpinOrder1,Sectral,Sectral1);
cout<<abs(trace(diagmat(RealizationError))/Realization);

Error<<x1<<"   "<<RealizationError(x1);
Error<<"   "<<abs(trace(diagmat(RealizationError))/(x1+1.00))<<endl;;


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

return 0;
}





