#include"Header.h"
void OptimizeIso(cx_mat*Oascend0,cx_mat *Uni0,cx_mat *Isometry0,cx_mat *Oascend1,cx_mat *Uni1,cx_mat *Isometry1){
cx_mat * Oout=new cx_mat[1];
cx_mat eigvec;
cx_vec eigval;
cx_mat U;
cx_mat V;
vec s;
double Gamma=1;
double E1,E2;
double E1h,E2h;
double NormReal;
cx_mat DX;
int Count;
int BreakWhile=1;
cx_mat Unih;
cx_mat IsoLh;
cx_mat Isoh;
cx_mat Z;
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
vec Xlvec2;
Xlvec2.zeros(6);
Xlvec2(0)=Xi2[0];
Xlvec2(1)=Xi2[1];
Xlvec2(2)=Xi2[2];
Xlvec2(3)=Xi2[3];
Xlvec2(4)=Xi2[4];
Xlvec2(5)=Xi2[5];

vec Xlvec3;
Xlvec3.zeros(6);
Xlvec3(0)=Xi3[0];
Xlvec3(1)=Xi3[1];
Xlvec3(2)=Xi3[2];
Xlvec3(3)=Xi3[3];
Xlvec3(4)=Xi3[4];
Xlvec3(5)=Xi3[5];


vec Xlvec4;
Xlvec4.zeros(6);
Xlvec4(0)=Xi4[0];
Xlvec4(1)=Xi4[1];
Xlvec4(2)=Xi4[2];
Xlvec4(3)=Xi4[3];
Xlvec4(4)=Xi4[4];
Xlvec4(5)=Xi4[5];


vec Xlvec5;
Xlvec5.zeros(6);
Xlvec5(0)=Xi5[0];
Xlvec5(1)=Xi5[1];
Xlvec5(2)=Xi5[2];
Xlvec5(3)=Xi5[3];
Xlvec5(4)=Xi5[4];
Xlvec5(5)=Xi5[5];

int NumIntChosen;
int NumInter=1;
int NumIso=0;
int NumUni=0;
int BoundryRight=1;
int BoundryLeft=1;
int NumInterLeft=0;
int NumInterRight=4;
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
cx_mat * Oascend1H=new cx_mat[5];
cx_mat * Oascend0H=new cx_mat[int(N0)];
cx_mat * Oascend2=new cx_mat[2];
cx_mat * Rho1=new cx_mat[int(N0/3)-1];
cx_mat * Rho2=new cx_mat[1];
cx_mat * Env=new cx_mat[1];
cx_mat * EnvU=new cx_mat[4];
cx_mat * IsoHelp=new cx_mat[int(N0/3)];
cx_mat * EnvIso=new cx_mat[2];
vec Xlvec0H=Xlvec0;
vec Xlvec2H=Xlvec2;


for(int j1=0;j1<1;j1++){



for(int q=2;q<=3;q++){
/********************************Isometry*****************************************/

if(q==3)
{
Xlvec2=Xlvec0H;
Xlvec0=Xlvec3;
}else{
Xlvec0=Xlvec0H;
Xlvec2=Xlvec2H;
}



NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec2,Oout);
Oascend1[q-1]=Oout[0];



NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=1;
BoundryLeft=1;
NumInterLeft=0;
NumInterRight=4;
for(int i=0;i<5;i++)
Oascend1H[i]=Oascend1[i];
Oascend1H[q-1].zeros(8*8,8*8);
Oascend1H[q].zeros(8*8,8*8);
Ascen(Oascend1H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
/*******************************/
NumIntChosen=q-1;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oout[0];
NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oascend2[1]+Oout[0];
/**************************************/

Oout[0]=Oascend2[1]+Oascend2[0];
cout<<"E"<<q<<"="<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;
Gamma=1;
Count=0;
for(int i=0;i<=NumIterationIso0;i++){

if(i==0)
Isoh=Isometry0[q];

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec0,Oout);
Oascend1[q]=Oout[0];


NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec2,Oout);
Oascend1[q-1]=Oout[0];


NumIntChosen=q-1;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oout[0];
NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oascend2[1]+Oout[0];

Oout[0]=Oascend2[1]+Oascend2[0];
//cout<<"E"<<q<<"="<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;

/*****************Energy**********************************/
Rho2[0]=Oascend2[1]+2*Oascend2[0];

NumRhoup=0;
NumRhodown=1;
NumRhoLeft=0;
NumRhoRight=4;
NumIntChosen=q-1;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumIntChosen=q;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

//cout<<"E"<<q<<"="<<((real(trace(Rho1[q-1]*Oascend1[q-1]))+real(trace(Rho1[q]*Oascend1[q]))
//+real((trace(Oascend2[0]*strans(Oascend2[0])))))*(1.00/pow(2,Nreal)))<<endl;

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumRhoup0=q;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoLeftMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryLeft0,
NumInterLeft0,Xlvec0,Oout);
EnvIso[0]=Oout[0];


NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
NumRhoup0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoRightMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryRight0,
NumInterRight0,Xlvec2,Oout);
EnvIso[1]=Oout[0];




//cout<<"E"<<q<<"="<<((real(trace((EnvIso[0]+EnvIso[1])*Isometry0[q]))))*(1.00/pow(2,Nreal))+
//real((trace(Oascend2[0]*strans(Oascend2[0]))))*(1.00/pow(2,Nreal))<<endl;
E1=((real(trace((EnvIso[0]+EnvIso[1])*Isometry0[q]))));
E1h=E1;
/*********************************Deritative-1*********************************************/
Rho2[0]=Oascend2[1];
NumIntChosen=q-1;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumIntChosen=q;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumRhoup0=q;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoLeftMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryLeft0,
NumInterLeft0,Xlvec0,Oout);
EnvIso[0]=Oout[0];


NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
NumRhoup0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoRightMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryRight0,
NumInterRight0,Xlvec2,Oout);
EnvIso[0]=Oout[0]+EnvIso[0];
/******************************************************************************/

/*********************************Deritative-2*********************************************/
Rho2[0]=2*Oascend2[0];
NumIntChosen=q-1;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumIntChosen=q;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumRhoup0=q;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoLeftMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryLeft0,
NumInterLeft0,Xlvec0,Oout);
EnvIso[1]=Oout[0];


NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
NumRhoup0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
EnvIsoRightMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,BoundryRight0,
NumInterRight0,Xlvec2,Oout);
EnvIso[1]=Oout[0]+EnvIso[1];
/******************************************************************************/


DX=trans(2*EnvIso[0]+EnvIso[1]);
Z=-(Isometry0[q]*trans(DX)*Isometry0[q]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-7)
break;
//cout<<"Norm="<<NormReal<<endl;
BreakWhile=1;
while(BreakWhile==1){
Count=Count+1;
IsoHelp[q]=Isometry0[q]+2.00*Gamma*Z;
svd_econ(U, s, V, IsoHelp[q],'b',"std");
IsoHelp[q]=U*trans(V);
/*****************Energy**********************************/
if( abs(Gamma) >1.0e+8 ){
cout<<"Break1, q"<<q<<" i="<<i<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}
IsoHelp[q+1]=Isometry0[q+1];
IsoHelp[q-1]=Isometry0[q-1];
NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,IsoHelp,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];
NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,IsoHelp,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec2,Oout);
Oascend1[q-1]=Oout[0];

NumIntChosen=q-1;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oout[0];
NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oascend2[1]+Oout[0];

Oout[0]=Oascend2[1]+Oascend2[0];
E2=real((trace(Oout[0]*strans(Oout[0]))))-real(trace(Oascend2[0]*strans(Oascend2[0])));
//cout<<"E1"<<q<<"="<<real((trace(Oout[0]*strans(Oout[0]))))*(1.00/pow(2,Nreal))<<endl;
/************************************************************************************/
if(E1-E2<=-NormReal*Gamma)
Gamma=2.00*Gamma;
else{
BreakWhile=0;
}
}
BreakWhile=1;

while(BreakWhile==1){

Count=Count+1;
IsoHelp[q]=Isometry0[q]+Gamma*Z;
svd_econ(U, s, V, IsoHelp[q],'b',"std");
IsoHelp[q]=U*trans(V);

/*****************Energy**********************************/
IsoHelp[q+1]=Isometry0[q+1];
IsoHelp[q-1]=Isometry0[q-1];
NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,IsoHelp,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];
NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,IsoHelp,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec2,Oout);
Oascend1[q-1]=Oout[0];
NumIntChosen=q-1;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,
Xlvec,Oout);
Oascend2[1]=Oout[0];
NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,Xlvec,Oout);
Oascend2[1]=Oascend2[1]+Oout[0];

Oout[0]=Oascend2[1]+Oascend2[0];
E2=real((trace(Oout[0]*strans(Oout[0]))))-real(trace(Oascend2[0]*strans(Oascend2[0])));
//cout<<"E2"<<q<<"="<<real((trace(Oout[0]*strans(Oout[0]))))*(1.00/pow(2,Nreal))<<endl;
/************************************************************************************/

if( abs((0.5)*NormReal*Gamma) <1.0e-12 ||  abs(E1-E2)<1.0e-12 ){
cout<<"Break, q"<<q<<" i="<<i<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}


if(E1-E2>=-(0.5)*NormReal*Gamma)
Gamma=(0.5)*Gamma;
else{
BreakWhile=0;
}
}

IsoHelp[q]=Isometry0[q]+Gamma*Z;
svd_econ(U, s, V, IsoHelp[q],'b',"std");
Isometry0[q]=U*trans(V);


if(E1h>E2h || i==0){
Isoh=Isometry0[q];
E2h=E1h;}
else{
Isometry0[q]=Isoh;
cout<<"NotOptimized="<<i<<"  E1h="<<E1h<<"  E2h"<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"Count1="<<Count<<endl;
break;
}


}


AscendTotal(Oascend0,Uni0,Isometry0,Oascend1);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<setprecision(16)<<"E"<<q<<"="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;


NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumInter0=3*q-2;
NumIso0=q-1;
NumUni0=q-1;
BoundryRight0=0;
BoundryLeft0=0;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,
Xlvec2,Oout);
Oascend1[q-1]=Oout[0];


}


}


delete [] Oout;
delete [] Oascend2;
delete [] Oascend1H;
delete [] Rho1;
delete [] Rho2;
delete [] Env;
delete [] EnvU;
delete [] IsoHelp;
delete [] EnvIso;
delete [] Oascend0H;


}
