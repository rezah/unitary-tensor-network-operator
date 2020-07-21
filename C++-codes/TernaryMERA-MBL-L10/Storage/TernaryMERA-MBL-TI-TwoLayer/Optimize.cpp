#include"Header.h"
void OptimizeU(cx_mat*Oascend0,cx_mat *Uni0,cx_mat *Isometry0,cx_mat *Oascend1,cx_mat *Uni1,cx_mat *Isometry1){


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
cx_mat IsoRh;
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


vec Xlvec0H=Xlvec0;

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
int NumRhoLeft=0;
int NumRhoRight=4;
int w;

int NumInter0;
int NumIso0;
int NumUni0;
int BoundryRight0;
int BoundryLeft0;
int NumInterLeft0=0;
int NumInterRight0=16;
int NumRhoup0;
int NumRhodown0;
int NumRhoLeft0;
int NumRhoRight0;

cx_mat * Oascend1H=new cx_mat[5];
cx_mat * Oascend2=new cx_mat[2];

cx_mat * Rho1=new cx_mat[int(N0/3)-1];
cx_mat * Rho2=new cx_mat[1];

cx_mat * Env=new cx_mat[1];
cx_mat * EnvU=new cx_mat[4];

cx_mat * UniHelp=new cx_mat[5];
for(int j1=0;j1<1;j1++){
for(int q=1;q<=3;q++){
/*********************************Uni-StepeseDescentOptimization*******************************/
if(q==2)
Xlvec0=Xlvec0H;
else if(q==1)
Xlvec0=Xlvec2;
else if(q==3)
Xlvec0=Xlvec3;



NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumInterLeft0=0;
NumInterRight0=16;
BoundryRight0=0;
BoundryLeft0=0;
if(q==4)
BoundryRight0=1;
if(q==0)
BoundryLeft0=1;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

cout<<setprecision(16);
NumInter=1;
 NumIso=0;
 NumUni=0;
 BoundryRight=1;
 BoundryLeft=1;
 NumInterLeft=0;
 NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<"E"<<q<<"="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;


for(int i=0;i<5;i++)
Oascend1H[i]=Oascend1[i];
Oascend1H[q].zeros(8*8,8*8);
Ascen(Oascend1H,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];


NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,Xlvec,Oout);
Oascend2[1]=Oout[0];

/*****************************Optimization*************************************/
Gamma=1;
Count=0;
for(int i=0;i<=NumIterationUni0;i++){

if(i==0)
Unih=Uni0[q];

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumInterLeft0=0;
NumInterRight0=16;
if(q==4)
BoundryRight0=1;
if(q==0)
BoundryLeft0=1;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,Xlvec,Oout);
Oascend2[1]=Oout[0];

/**************Env********************/
Rho2[0]=2*Oascend2[0];
NumRhoup=0;
NumRhodown=1;
NumRhoLeft=0;
NumRhoRight=4;
NumIntChosen=q; /*****0-4******/
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumRhoup0=q;
EnvUniMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,Xlvec0,Oout);
EnvU[1]=Oout[0];
/*******************************/
Rho2[0]=Oascend2[1];
NumRhoup=0;
NumRhodown=1;
NumRhoLeft=0;
NumRhoRight=4;
NumIntChosen=q;
Dscenone(Rho2,Rho1, Uni1, Isometry1,NumRhoup,NumRhodown,NumIso,NumUni,NumRhoLeft,
NumRhoRight,NumIntChosen,Xlvec);
NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumRhoup0=q;
EnvUniMERA(Oascend0,Rho1, Uni0, Isometry0,NumInter0,NumRhoup0,NumIso0,NumUni0,Xlvec0,Oout);
EnvU[2]=Oout[0];
/**********************Energy*****************************/
Oout[0]=Oascend2[1]+Oascend2[0];
E1=real((trace(Oout[0]*strans(Oout[0]))))-real(trace(Oascend2[0]*strans(Oascend2[0])));
E1h=E1;
/**********************End-Energy*****************************/
/*********************Optimization*************************/
DX=trans(EnvU[1]+2*EnvU[2]);
Z=-(Uni0[q]*trans(DX)*Uni0[q]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-6){
cout<<"BreakNorm="<<NormReal<<endl;
break;}
//cout<<"Norm="<<NormReal<<endl;
BreakWhile=1;
while(BreakWhile==1){
Count=Count+1;
UniHelp[q]=Uni0[q]+2.00*Gamma*Z;
svd_econ(U, s, V, UniHelp[q],'b',"std");
UniHelp[q]=U*trans(V);

/**********************Energy*****************************/

if( abs(Gamma) >1.0e+8 ){
cout<<"Break1, q"<<q<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}

NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumInterLeft0=0;
NumInterRight0=16;
if(q==4)
BoundryRight0=1;
if(q==0)
BoundryLeft0=1;
AscenMERA(Oascend0,UniHelp,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,Xlvec,Oout);
Oascend2[1]=Oout[0];


Oout[0]=Oascend2[1]+Oascend2[0];
E2=real((trace(Oout[0]*strans(Oout[0]))))-real(trace(Oascend2[0]*strans(Oascend2[0])));
//cout<<"E2="<<E2<<endl;
//cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
//cout<<"Count1="<<Count<<endl;
/**********************Energy*****************************/
if(E1-E2<=-NormReal*Gamma)
Gamma=2.00*Gamma;
else{
BreakWhile=0;
}
}

BreakWhile=1;
while(BreakWhile==1){
Count=Count+1;
UniHelp[q]=Uni0[q]+Gamma*Z;
svd_econ(U, s, V, UniHelp[q],'b',"std");
UniHelp[q]=U*trans(V);


/**********************Energy*****************************/
NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumInterLeft0=0;
NumInterRight0=16;
if(q==4)
BoundryRight0=1;
if(q==0)
BoundryLeft0=1;
AscenMERA(Oascend0,UniHelp,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumIntChosen=q;
AscenOne(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,NumIntChosen,Xlvec,Oout);
Oascend2[1]=Oout[0];


Oout[0]=Oascend2[1]+Oascend2[0];
E2=real((trace(Oout[0]*strans(Oout[0]))))-real(trace(Oascend2[0]*strans(Oascend2[0])));


/**********************Energy*****************************/

if( abs((0.5)*NormReal*Gamma) <1.0e-12 ||  abs(E1-E2)<1.0e-12 ){
cout<<"Break2, q"<<q<<" i="<<i<<endl;
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

UniHelp[q]=Uni0[q]+Gamma*Z;
svd_econ(U, s, V, UniHelp[q],'b',"std");
Uni0[q]=U*trans(V);


if(E1h>E2h || i==0){
Unih=Uni0[q];
E2h=E1h;}
else{
Uni0[q]=Unih;
cout<<"NotOptimized, i="<<i<<" q="<<q<<"  E1h="<<E1h<<"  E2h"<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"q="<<q<<" Count1="<<Count<<endl;
break;
}

}


NumInter0=3*q+1;
NumIso0=q;
NumUni0=q;
NumInterLeft0=0;
NumInterRight0=16;
if(q==4)
BoundryRight0=1;
if(q==0)
BoundryLeft0=1;
AscenMERA(Oascend0,Uni0,Isometry0,NumInter0,NumIso0,NumUni0,BoundryRight0,BoundryLeft0,NumInterLeft0,NumInterRight0,Xlvec0,Oout);
Oascend1[q]=Oout[0];

NumInter=1;
NumIso=0;
NumUni=0;
BoundryRight=1;
BoundryLeft=1;
NumInterLeft=0;
NumInterRight=4;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
cout<<"E"<<q<<"="<<abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;

/*****************************End*************************************/
}
}
delete [] Oout;
delete [] Oascend2;
delete [] Oascend1H;
delete [] Rho1;
delete [] Rho2;
delete [] Env;
delete [] EnvU;
delete [] UniHelp;
}
