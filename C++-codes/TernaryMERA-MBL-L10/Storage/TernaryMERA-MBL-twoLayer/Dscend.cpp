#include"Header.h"
void Dscen(cx_mat*Rho2,cx_mat*Rho1,cx_mat*Uni,cx_mat*Isometry,
int NumRhoup,int NumRhodown,int NumIso,int NumUni,int BoundryRight,int BoundryLeft,
int NumRhoLeft,int NumRhoRight,vec Xlvec)
{
    complex<double> sum;


int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);

    cx_mat Ohelp;
    cx_mat Ohelp1;
    cx_mat Ohelp2;
/*************************B-Midle*****************************/
Ohelp.zeros(X2*X2,X0*X1*X2);

for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp1=Ohelp*Rho2[NumRhoup];

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp1=Ohelp1*Ohelp;


Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp1(mb*X2+m,xb*X2b+x);

Ohelp=Uni[NumUni]*Ohelp*trans(Uni[NumUni]);


Rho1[NumRhodown+1]=Ohelp;
/***************************B-Left********************************************/

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp1=Rho2[NumRhoup]*Ohelp;

Ohelp.zeros(X2*X2*X1*X1,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,q));
Ohelp(jb*X1*X2*X2+j*X2*X2+mb*X2+m,q)=sum;
}
Ohelp2=Ohelp*Ohelp1;



Ohelp1.zeros(X2b*X2*X2*X2b,X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j)=
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x);


Ohelp.zeros(X2b,X2b*X2*X2);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2*X2b+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;
Ohelp2.zeros(X2*X2,X2b*X2b*X2*X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(sb*X2+s,xb*X2b*X2*X2+x*X2*X2+mb*X2+m)=
Ohelp(sb*X2*X2b+mb*X2b+xb,s*X2b*X2+m*X2b+x);

Ohelp1=Ohelp2*Ohelp1;

Ohelp.zeros(X2*X1,X2*X1);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(j*X2+s,jb*X2+sb)=Ohelp1(sb*X2+s,jb*X1+j);

Rho1[NumRhodown]=Ohelp;
////cout<<endl<<setprecision(16)<<trace(Ohelp*Oascend[NumLay-1])<<endl;
/*******************B-Right************************/
Ohelp.zeros(X2*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp1=Ohelp*Rho2[NumRhoup];

Ohelp.zeros(X3*X4*X2b,X2b*X2b*X3*X3);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,g));
Ohelp(g,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x)=sum;
}
Ohelp2=Ohelp1*Ohelp;


Ohelp1.zeros(X2b*X2*X2*X2b,X3*X3);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v)=
Ohelp2(mb*X2+m,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x);


Ohelp.zeros(X2,X2b*X2b*X2);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(s,c*X2*X2b+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;

Ohelp2.zeros(X2b*X2b,X2b*X2b*X2*X2);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(cb*X2b+c,xb*X2b*X2*X2+x*X2*X2+mb*X2+m)=
Ohelp(cb*X2*X2b+mb*X2b+xb,c*X2b*X2+m*X2b+x);

Ohelp1=Ohelp2*Ohelp1;
Ohelp.zeros(X2b*X3,X2b*X3);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(c*X3+v,cb*X3+vb)=Ohelp1(cb*X2b+c,vb*X3+v);

Rho1[NumRhodown+2]=Ohelp;
/***************************BoundryLeft***********************************/

if(BoundryLeft==1){

Ohelp1.zeros(X0*X1*X2,1);
for(int q=0;q<X0*X1*X2;q++){sum=0;
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Rho2[NumRhoup](q,g);
Ohelp1(q,0)=sum;
}

Ohelp.zeros(X0*X1,X0*X1);
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int m=0;m<X2;m++)
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Ohelp1(q,0)*conj(Isometry[NumIso](ib*X2*X1+jb*X2+m,q));
Ohelp(i*X1+j,ib*X1+jb)=sum;
}
Rho1[NumRhoLeft]=Ohelp;
}

/***************************Boundryight****************************************/
if(BoundryRight==1){

Ohelp1.zeros(X3*X4*X2b,1);
for(int g=0;g<X3*X4*X2b;g++){sum=0;
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Rho2[NumRhoup](q,g);
Ohelp1(g,0)=sum;
}

Ohelp.zeros(X3*X4,X3*X4);
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int x=0;x<X2b;x++)
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Ohelp1(g,0)*conj(Isometry[NumIso+1](x*X3*X4+vb*X4+zb,g));
Ohelp(v*X4+z,vb*X4+zb)=sum;
}
Rho1[NumRhoRight]=Ohelp;
}



}


void Dscenone(cx_mat*Rho2,cx_mat*Rho1,cx_mat*Uni,cx_mat*Isometry,
int NumRhoup,int NumRhodown,int NumIso,int NumUni,
int NumRhoLeft,int NumRhoRight, int IntnumChosen ,vec Xlvec)
{
    complex<double> sum;


int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);

    cx_mat Ohelp;
    cx_mat Ohelp1;
    cx_mat Ohelp2;
if(IntnumChosen==2){
/*************************B-Midle*****************************/
Ohelp.zeros(X2*X2,X0*X1*X2);

for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp1=Ohelp*Rho2[NumRhoup];

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp1=Ohelp1*Ohelp;


Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp1(mb*X2+m,xb*X2b+x);

Ohelp=Uni[NumUni]*Ohelp*trans(Uni[NumUni]);


Rho1[NumRhodown+1]=Ohelp;}

if(IntnumChosen==1){
/***************************B-Left********************************************/

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp1=Rho2[NumRhoup]*Ohelp;

Ohelp.zeros(X2*X2*X1*X1,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,q));
Ohelp(jb*X1*X2*X2+j*X2*X2+mb*X2+m,q)=sum;
}
Ohelp2=Ohelp*Ohelp1;



Ohelp1.zeros(X2b*X2*X2*X2b,X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j)=
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x);


Ohelp.zeros(X2b,X2b*X2*X2);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2*X2b+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;
Ohelp2.zeros(X2*X2,X2b*X2b*X2*X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(sb*X2+s,xb*X2b*X2*X2+x*X2*X2+mb*X2+m)=
Ohelp(sb*X2*X2b+mb*X2b+xb,s*X2b*X2+m*X2b+x);

Ohelp1=Ohelp2*Ohelp1;

Ohelp.zeros(X2*X1,X2*X1);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(j*X2+s,jb*X2+sb)=Ohelp1(sb*X2+s,jb*X1+j);

Rho1[NumRhodown]=Ohelp;}
////cout<<endl<<setprecision(16)<<trace(Ohelp*Oascend[NumLay-1])<<endl;

if(IntnumChosen==3){
/*******************B-Right************************/
Ohelp.zeros(X2*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp1=Ohelp*Rho2[NumRhoup];

Ohelp.zeros(X3*X4*X2b,X2b*X2b*X3*X3);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,g));
Ohelp(g,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x)=sum;
}
Ohelp2=Ohelp1*Ohelp;


Ohelp1.zeros(X2b*X2*X2*X2b,X3*X3);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v)=
Ohelp2(mb*X2+m,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x);


Ohelp.zeros(X2,X2b*X2b*X2);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(s,c*X2*X2b+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;

Ohelp2.zeros(X2b*X2b,X2b*X2b*X2*X2);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(cb*X2b+c,xb*X2b*X2*X2+x*X2*X2+mb*X2+m)=
Ohelp(cb*X2*X2b+mb*X2b+xb,c*X2b*X2+m*X2b+x);

Ohelp1=Ohelp2*Ohelp1;
Ohelp.zeros(X2b*X3,X2b*X3);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(c*X3+v,cb*X3+vb)=Ohelp1(cb*X2b+c,vb*X3+v);

Rho1[NumRhodown+2]=Ohelp;}
/***************************BoundryLeft***********************************/

if(IntnumChosen==0){

Ohelp1.zeros(X0*X1*X2,1);
for(int q=0;q<X0*X1*X2;q++){sum=0;
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Rho2[NumRhoup](q,g);
Ohelp1(q,0)=sum;
}

Ohelp.zeros(X0*X1,X0*X1);
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int m=0;m<X2;m++)
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Ohelp1(q,0)*conj(Isometry[NumIso](ib*X2*X1+jb*X2+m,q));
Ohelp(i*X1+j,ib*X1+jb)=sum;
}
Rho1[NumRhoLeft]=Ohelp;
}

/***************************Boundryight****************************************/
if( IntnumChosen==4){

Ohelp1.zeros(X3*X4*X2b,1);
for(int g=0;g<X3*X4*X2b;g++){sum=0;
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Rho2[NumRhoup](q,g);
Ohelp1(g,0)=sum;
}

Ohelp.zeros(X3*X4,X3*X4);
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int x=0;x<X2b;x++)
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Ohelp1(g,0)*conj(Isometry[NumIso+1](x*X3*X4+vb*X4+zb,g));
Ohelp(v*X4+z,vb*X4+zb)=sum;
}
Rho1[NumRhoRight]=Ohelp;
}



}


