#include"Header.h"
void OptimizeUp(cx_mat*Oascend0,cx_mat *Uni0,cx_mat *Isometry0,
cx_mat *Oascend1,cx_mat *Uni1,cx_mat *Isometry1,vec *Checker,cx_mat *Oascend2){

cx_mat U;
cx_mat V;
vec s;

int NumInter;
int NumIso;
int NumUni;vec Xlvec;

int BoundryRight;
int BoundryLeft;
int NumInterLeft;
int NumInterRight;
int NumRhoup;
int NumRhodown;
int NumRhoLeft;
int NumRhoRight;
int w;

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
 NumRhoup=0;
 NumRhodown=1;
 NumRhoLeft=0;
 NumRhoRight=4;
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

cx_mat * Rho2=new cx_mat[1];
cx_mat * EnvU=new cx_mat[4];
cx_mat * Oout=new cx_mat[1];
cx_mat * EnvIso=new cx_mat[4];

/***********************************Uni**************************************************/
Gamma=1;
for(int i=0;i<=NumIterationUP;i++){
if(i==0)
Unih=Uni1[0];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E1=real(trace(Oascend2[0]*strans(Oascend2[0])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));
E1h=E1;
//cout<<"E="<<(E1+real(trace(Oascend2[1]*strans(Oascend2[1]))))*(1.00/pow(2,Nreal))<<endl;
Rho2[0]=Oascend2[1]+Oascend2[0];
NumInter=1;
NumIso=0;
NumUni=0;
EnvUni(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
Xlvec,Oout);
DX=trans(Oout[0]);
Z=-(Uni1[0]*trans(DX)*Uni1[0]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-7)
break;
BreakWhile=1;
while(BreakWhile==1){
if( abs(Gamma) >1.0e+8 ){
cout<<"Break1"<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}
Count=Count+1;
EnvU[0]=Uni1[0]+2.00*Gamma*Z;
svd_econ(U, s, V, EnvU[0],'b',"std");
EnvU[0]=U*trans(V);
EnvU[1]=Uni1[1];

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E2=real(trace(Oascend2[0]*strans(Oascend2[0])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));
if(E1-E2<=-NormReal*Gamma)
Gamma=2.00*Gamma;
else{
BreakWhile=0;
}


}
BreakWhile=1;


while(BreakWhile==1){
Count=Count+1;
EnvU[0]=Uni1[0]+Gamma*Z;
svd_econ(U, s, V, EnvU[0],'b',"std");
EnvU[0]=U*trans(V);
EnvU[1]=Uni1[1];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E2=real(trace(Oascend2[0]*strans(Oascend2[0])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));



if( abs((0.5)*NormReal*Gamma) <1.0e-13 ||  abs(E1-E2)<1.0e-12 ){
i=i+120;
cout<<"Break"<<endl;
break;
}


if(E1-E2>=-(0.5)*NormReal*Gamma)
Gamma=(0.5)*Gamma;
else{
BreakWhile=0;
}
}


EnvU[0]=Uni1[0]+Gamma*Z;
svd_econ(U, s, V, EnvU[0],'b',"std");
Uni1[0]=U*trans(V);

if(E1h>E2h || i==0){
Unih=Uni1[0];
E2h=E1h;}
else{
Uni1[0]=Unih;
cout<<"NotOptimizedUniup="<<i<<"  E1h="<<E1h<<"  E2h"<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"Count1="<<Count<<endl;
break;
}


}

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oout[0]=strans(Oout[0])+Oascend2[0];
cout<<"E1="<<setprecision(16)<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;



/***********************************Uni2**************************************************/
Gamma=1;
for(int i=0;i<=NumIterationUP;i++){
if(i==0)
Unih=Uni1[1];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E1=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));
E1h=E1;
//cout<<"E="<<(E1+real(trace(Oascend2[0]*strans(Oascend2[0]))))*(1.00/pow(2,Nreal))<<endl;
Rho2[0]=strans(Oascend2[1])+strans(Oascend2[0]);
NumInter=4;
NumIso=1;
NumUni=1;
EnvUni(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,
Xlvec,Oout);
DX=trans(Oout[0]);
Z=-(Uni1[1]*trans(DX)*Uni1[1]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-7)
break;
BreakWhile=1;
while(BreakWhile==1){
if( abs(Gamma) >1.0e+8 ){
cout<<"Break1"<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}
Count=Count+1;
EnvU[1]=Uni1[1]+2.00*Gamma*Z;
svd_econ(U, s, V, EnvU[1],'b',"std");
EnvU[1]=U*trans(V);
EnvU[0]=Uni1[0];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));

if(E1-E2<=-NormReal*Gamma)
Gamma=2.00*Gamma;
else{
BreakWhile=0;
}


}
BreakWhile=1;

while(BreakWhile==1){
Count=Count+1;
EnvU[1]=Uni1[1]+Gamma*Z;
svd_econ(U, s, V, EnvU[1],'b',"std");
EnvU[1]=U*trans(V);
EnvU[0]=Uni1[0];

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,EnvU,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])));

if( abs((0.5)*NormReal*Gamma) <1.0e-13 ||  abs(E1-E2)<1.0e-12 ){
i=i+120;
cout<<"Break"<<endl;
break;
}


if(E1-E2>=-(0.5)*NormReal*Gamma)
Gamma=(0.5)*Gamma;
else{
BreakWhile=0;
}
}


EnvU[1]=Uni1[1]+Gamma*Z;
svd_econ(U, s, V, EnvU[1],'b',"std");
Uni1[1]=U*trans(V);

if(E1h>E2h || i==0){
Unih=Uni1[1];
E2h=E1h;}
else{
Uni1[1]=Unih;
cout<<"NotOptimizedUniup="<<i<<"  E1h="<<E1h<<"  E2h="<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"Count1="<<Count<<endl;
break;
}


}
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oout[0]=strans(Oout[0])+Oascend2[0];
cout<<"E2="<<setprecision(16)<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;

/*************************************Left*****************************************/
Gamma=1;
for(int i=0;i<=NumIterationUP;i++){
if(i==0)
Unih=Isometry1[0];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E1=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));
E1h=E1;
//cout<<"E="<<(E1)*(1.00/pow(2,Nreal))<<endl;

NumInter=1;
NumIso=0;
NumUni=0;
Rho2[0]=Oascend2[0]+Oascend2[1];
EnvIsoLeft(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,
NumInterLeft,NumInterRight,Xlvec,Oout);
EnvU[0]=Oout[0];

Rho2[0]=strans(Oascend2[1])+strans(Oascend2[0]);
NumInter=4;
NumIso=1;
NumUni=1;
//cout<<BoundryRight<<BoundryRight<<endl;
EnvIsoRight(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,
NumInterLeft,NumInterRight,Xlvec,Oout);
EnvU[1]=Oout[0];
//cout<<"ET="<<real(trace((EnvU[1]+EnvU[0])*(Isometry1[0])))*(1.00/pow(2,Nreal))<<endl;

DX=trans(EnvU[1]+EnvU[0]);
Z=-(Isometry1[0]*trans(DX)*Isometry1[0]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-8)
break;
//cout<<endl<<"NormReal="<<NormReal<<endl;
Count=0;
BreakWhile=1;
while(BreakWhile==1){

if( abs(Gamma) >1.0e+8 ){
cout<<"Break1"<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}
Count=Count+1;
EnvIso[0]=Isometry1[0]+2*Gamma*Z;
svd_econ(U, s, V, EnvIso[0],'b',"std");
EnvIso[0]=U*trans(V);
EnvIso[1]=Isometry1[1];
EnvIso[2]=EnvIso[0];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);

E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));
//cout<<"E21="<<(E2)*(1.00/pow(2,Nreal))<<endl;


if(E1-E2<=-NormReal*Gamma)
Gamma=2*Gamma;
else{
//cout<<E1-E2<<" 1 "<< NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
//cout<<"Count1="<<Count<<endl;
BreakWhile=0;
}
}

BreakWhile=1;
Count=0;
while(BreakWhile==1){
Count=Count+1;
EnvIso[0]=Isometry1[0]+Gamma*Z;
svd_econ(U, s, V, EnvIso[0],'b',"std");
EnvIso[0]=U*trans(V);
EnvIso[1]=Isometry1[1];
EnvIso[2]=EnvIso[0];

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);

E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));


if( abs((0.5)*NormReal*Gamma) <1.0e-12 ||  abs(E1-E2)<1.0e-12 ){
i=i+120;
cout<<"Break"<<endl;
cout<<"Break"<<" i="<<i<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
break;
}

if(E1-E2>=-(0.5)*NormReal*Gamma)
Gamma=(0.5)*Gamma;
else{
BreakWhile=0;
//cout<<E1-E2<<" 2 "<< -(0.5)*NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
//cout<<"Count2="<<Count<<endl;
}
}
EnvIso[0]=Isometry1[0]+Gamma*Z;
svd_econ(U, s, V, EnvIso[0],'b',"std");
Isometry1[0]=U*trans(V);
Isometry1[2]=Isometry1[0];

if(E1h>E2h || i==0){
Unih=Isometry1[0];
Isometry1[2]=Isometry1[0];
E2h=E1h;}
else{
Uni1[0]=Isometry1[0];
Isometry1[2]=Isometry1[0];
cout<<"NotOptimizedIso1up="<<i<<"  E1h="<<E1h<<"  E2h"<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"Count1="<<Count<<endl;
break;
}

}

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oout[0]=strans(Oout[0])+Oascend2[0];
cout<<"E3="<<setprecision(16)<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;
/************************************Right**********************************************************************/
Gamma=1;
for(int i=0;i<=NumIterationUP;i++){
if(i==0)
Unih=Isometry1[1];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);
E1=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));
E1h=E1;
//cout<<"E="<<(E1)*(1.00/pow(2,Nreal))<<endl;

NumInter=1;
NumIso=0;
NumUni=0;
Rho2[0]=Oascend2[0]+Oascend2[1];
EnvIsoRight(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,
NumInterLeft,NumInterRight,Xlvec,Oout);
EnvU[0]=Oout[0];

Rho2[0]=strans(Oascend2[1])+strans(Oascend2[0]);
NumInter=4;
NumIso=1;
NumUni=1;
//cout<<BoundryRight<<BoundryRight<<endl;
EnvIsoLeft(Oascend1,Rho2, Uni1, Isometry1,NumInter,NumRhoup,NumIso,NumUni,BoundryRight,BoundryLeft,
NumInterLeft,NumInterRight,Xlvec,Oout);
EnvU[1]=Oout[0];
//cout<<"ET="<<real(trace((EnvU[1]+EnvU[0])*(Isometry1[0])))*(1.00/pow(2,Nreal))<<endl;

DX=trans(EnvU[1]+EnvU[0]);
Z=-(Isometry1[1]*trans(DX)*Isometry1[1]-DX);
NormReal=real(trace(trans(Z)*Z)*0.5);
if(NormReal<1.0e-8)
break;
//cout<<endl<<"NormReal="<<NormReal<<endl;
Count=0;
BreakWhile=1;
while(BreakWhile==1){

if( abs(Gamma) >1.0e+8 ){
cout<<"Break1"<<endl;
cout<<"NormReal="<<NormReal<<"  Gamma= "<<Gamma<<endl;
cout<<"Count="<<Count<<endl;
i=i+100;
break;
}
Count=Count+1;
EnvIso[1]=Isometry1[1]+2*Gamma*Z;
svd_econ(U, s, V, EnvIso[1],'b',"std");
EnvIso[1]=U*trans(V);
EnvIso[0]=Isometry1[0];
EnvIso[2]=EnvIso[0];
NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);

E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));
//cout<<"E21="<<(E2)*(1.00/pow(2,Nreal))<<endl;


if(E1-E2<=-NormReal*Gamma)
Gamma=2*Gamma;
else{
//cout<<E1-E2<<" 1 "<< NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
//cout<<"Count1="<<Count<<endl;
BreakWhile=0;
}
}

BreakWhile=1;
Count=0;
while(BreakWhile==1){
Count=Count+1;
EnvIso[1]=Isometry1[1]+Gamma*Z;
svd_econ(U, s, V, EnvIso[1],'b',"std");
EnvIso[1]=U*trans(V);
EnvIso[0]=Isometry1[0];
EnvIso[2]=EnvIso[0];

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,EnvIso,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[1]=strans(Oout[0]);

E2=real(trace(Oascend2[1]*strans(Oascend2[1])))+2*real(trace(Oascend2[0]*strans(Oascend2[1])))+
real(trace(Oascend2[0]*strans(Oascend2[0])));


if( abs((0.5)*NormReal*Gamma) <1.0e-12 ||  abs(E1-E2)<1.0e-12 ){
i=i+120;
cout<<"Break"<<endl;
break;
}

if(E1-E2>=-(0.5)*NormReal*Gamma)
Gamma=(0.5)*Gamma;
else{
BreakWhile=0;
//cout<<E1-E2<<" 2 "<< -(0.5)*NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
//cout<<"Count2="<<Count<<endl;
}
}
EnvIso[1]=Isometry1[1]+Gamma*Z;
svd_econ(U, s, V, EnvIso[1],'b',"std");
Isometry1[1]=U*trans(V);

if(E1h>E2h || i==0){
Unih=Isometry1[1];
E2h=E1h;}
else{
Uni1[0]=Isometry1[1];
cout<<"NotOptimizedIso1up="<<i<<"  E1h="<<E1h<<"  E2h"<<E2h<<endl;
cout<<E1-E2<<" 1 "<< -NormReal*Gamma<<"  Gamma= "<<Gamma<<endl;
cout<<"Count1="<<Count<<endl;
break;
}

}

NumInter=1;
NumIso=0;
NumUni=0;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oascend2[0]=Oout[0];
NumInter=4;
NumIso=1;
NumUni=1;
Ascen(Oascend1,Uni1,Isometry1,NumInter,NumIso,NumUni,BoundryRight,BoundryLeft,NumInterLeft,NumInterRight,Xlvec,Oout);
Oout[0]=strans(Oout[0])+Oascend2[0];
cout<<"E4="<<setprecision(16)<<real((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)))<<endl;

Checker[0](0)=Checker[0](1);
Checker[0](1)=abs((trace(Oout[0]*strans(Oout[0])))*(1.00/pow(2,Nreal)));



delete [] Rho2;
delete [] EnvU;
delete [] Oout;
delete [] EnvIso;

}
