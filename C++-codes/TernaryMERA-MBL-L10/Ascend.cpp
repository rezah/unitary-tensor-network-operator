#include"Header.h"
void Ascen(cx_mat*Oascend,cx_mat*Uni,cx_mat*Isometry,
int NumInter,int NumIso,int NumUni,int BoundryRight,int BoundryLeft,int NumInterLeft,int NumInterRight,
vec Xlvec,cx_mat *Oout){
cx_mat Ohelp;
    cx_mat Ohelp1;
    cx_mat Ohelp2;
    cx_mat Ohelp3;
int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);
complex<double> sum;
/*************************B-Midle*****************************/
Ohelp.zeros(X0*X1*X2,X2*X2);


for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(q,mb*X2+m)=sum;
}




Ohelp2=trans(Uni[NumUni])*Oascend[NumInter+1]*Uni[NumUni];



Ohelp1.zeros(X2*X2,X2b*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp1(mb*X2+m,xb*X2b+x)=
Ohelp2(mb*X2b+xb,m*X2b+x);

Ohelp1=Ohelp*Ohelp1;

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Ohelp1=Ohelp1*strans(Ohelp);

Oout[0]=Ohelp1;
//cout<<setprecision(16)<<trace(Oout[0]*strans(Oout[0]))<<endl;
/************************B-Left**************************/
Ohelp.zeros(X2b,X2*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);


Ohelp=trans(Ohelp)*Ohelp;

Ohelp1.zeros(X2b*X2b*X2*X2,X2*X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,sb*X2+s)=Ohelp(sb*X2b*X2+mb*X2b+xb,s*X2b*X2+m*X2b+x);


Ohelp.zeros(X2*X2,X1*X1);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(sb*X2+s,jb*X1+j)=Oascend[NumInter](jb*X2+sb,j*X2+s);

Ohelp=Ohelp1*Ohelp;
Ohelp2.zeros(X2*X2*X1*X1,X2b*X2b);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j);

Ohelp.zeros(X2*X1*X0,X2*X2*X1*X1);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,q));
Ohelp(q,jb*X1*X2*X2+j*X2*X2+mb*X2+m)=sum;
}
Ohelp2=Ohelp*Ohelp2;


Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Oout[0]=Oout[0]+Ohelp2*strans(Ohelp);
//Oout[0]=Ohelp2*strans(Ohelp);
//cout<<setprecision(16)<<trace(Oout[0]*strans(Oout[0]))<<endl;

/************************B-Right**************************/
Ohelp.zeros(X2,X2b*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(s,c*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;


Ohelp1.zeros(X2b*X2b*X2*X2,X2b*X2b);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,cb*X2b+c)=Ohelp(cb*X2b*X2+mb*X2b+xb,c*X2b*X2+m*X2b+x);

Ohelp.zeros(X2b*X2b,X3*X3);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(cb*X2b+c,vb*X3+v)=Oascend[NumInter+2](cb*X3+vb,c*X3+v);


Ohelp=Ohelp1*Ohelp;
Ohelp2.zeros(X2b*X2b*X3*X3,X2*X2);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v);


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

Ohelp2=strans(Ohelp2)*strans(Ohelp);

Ohelp.zeros(X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(q,mb*X2+m)=sum;
}
Oout[0]=Oout[0]+Ohelp*Ohelp2;


/*****************************Left******************************/
if(BoundryLeft==1){

Ohelp.zeros(X0*X1*X2,X0*X1*X2);


for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Oascend[NumInterLeft](ib*X1+jb,i*X1+j);
Ohelp(ib*X2*X1+jb*X2+m,q)=sum;
}


Ohelp1.zeros(X0*X1*X2,X2b*X3*X4);
for(int q=0;q<X0*X1*X2;q++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int m=0;m<X2;m++)
sum=sum+Ohelp(i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+m,q));
for(int qb=0;qb<X2b*X3*X4;qb++)
Ohelp1(q,qb)=sum;
}

Oout[0]=Oout[0]+Ohelp1;

}
/*****************************Right******************************/

if(BoundryRight==1){

Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Oascend[NumInterRight](vb*X4+zb,v*X4+z);
Ohelp(x*X3*X4+vb*X4+zb,g)=sum;
}


Ohelp1.zeros(X0*X1*X2,X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
for(int x=0;x<X2b;x++)
sum=sum+Ohelp(x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](x*X3*X4+v*X4+z,g));
for(int qb=0;qb<X2*X1*X0;qb++)
Ohelp1(qb,g)=sum;
}

Oout[0]=Oout[0]+Ohelp1;
//cout<<setprecision(16)<<trace(Oout[0]*strans(Oout[0]))<<endl;


}


}


void AscenMERA(cx_mat*Oascend,cx_mat*Uni,cx_mat*Isometry,
int NumInter,int NumIso,int NumUni,int BoundryRight,int BoundryLeft,int NumInterLeft,int NumInterRight,
vec Xlvec,cx_mat *Oout){
cx_mat Ohelp;
    cx_mat Ohelp1;
    cx_mat Ohelp2;
    cx_mat Ohelp3;
int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);
complex<double> sum;
/*************************B-Midle*****************************/
Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);

for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(q*X0*X1*X2+qb,mb*X2+m)=sum;
}



Ohelp2=trans(Uni[NumUni])*Oascend[NumInter+1]*Uni[NumUni];


Ohelp1.zeros(X2*X2,X2b*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp1(mb*X2+m,xb*X2b+x)=
Ohelp2(mb*X2b+xb,m*X2b+x);

Ohelp1=Ohelp*Ohelp1;

Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(g*X3*X4*X2b+gb,xb*X2b+x)=sum;
}
Ohelp1=Ohelp1*strans(Ohelp);


Ohelp.zeros(X2*X0*X1*X2b*X3*X4,X2*X0*X1*X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp(qb*X3*X4*X2b+gb,q*X3*X4*X2b+g)=Ohelp1(q*X0*X1*X2+qb,g*X3*X4*X2b+gb);


Oout[0]=Ohelp;
/************************B-Left**************************/
Ohelp.zeros(X2b,X2*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);
Ohelp=trans(Ohelp)*Ohelp;

Ohelp1.zeros(X2b*X2b*X2*X2,X2*X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,sb*X2+s)=Ohelp(sb*X2b*X2+mb*X2b+xb,s*X2b*X2+m*X2b+x);
Ohelp.zeros(X2*X2,X1*X1);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(sb*X2+s,jb*X1+j)=Oascend[NumInter](jb*X2+sb,j*X2+s);

Ohelp=Ohelp1*Ohelp;
Ohelp2.zeros(X2*X2*X1*X1,X2b*X2b);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j);

Ohelp.zeros(X2*X1*X0*X2*X1*X0,X2*X2*X1*X1);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,qb));
Ohelp(q*X2*X1*X0+qb,jb*X1*X2*X2+j*X2*X2+mb*X2+m)=sum;
}
Ohelp2=Ohelp*Ohelp2;


Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(g*X3*X4*X2b+gb,xb*X2b+x)=sum;
}


Ohelp2=Ohelp2*strans(Ohelp);

Ohelp.zeros(X2*X0*X1*X2b*X3*X4,X2*X0*X1*X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp(qb*X3*X4*X2b+gb,q*X3*X4*X2b+g)=Ohelp2(q*X0*X1*X2+qb,g*X3*X4*X2b+gb);


Oout[0]=Oout[0]+Ohelp;
/************************B-Right**************************/
Ohelp.zeros(X2,X2b*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(s,c*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;


Ohelp1.zeros(X2b*X2b*X2*X2,X2b*X2b);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,cb*X2b+c)=Ohelp(cb*X2b*X2+mb*X2b+xb,c*X2b*X2+m*X2b+x);

Ohelp.zeros(X2b*X2b,X3*X3);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(cb*X2b+c,vb*X3+v)=Oascend[NumInter+2](cb*X3+vb,c*X3+v);


Ohelp=Ohelp1*Ohelp;
Ohelp2.zeros(X2b*X2b*X3*X3,X2*X2);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v);


Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b*X3*X3);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,gb));
Ohelp(g*X3*X4*X2b+gb,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x)=sum;
}

Ohelp2=strans(Ohelp2)*strans(Ohelp);

Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(q*X0*X1*X2+qb,mb*X2+m)=sum;
}

Ohelp2=Ohelp*Ohelp2;

Ohelp.zeros(X2*X0*X1*X2b*X3*X4,X2*X0*X1*X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp(qb*X3*X4*X2b+gb,q*X3*X4*X2b+g)=Ohelp2(q*X0*X1*X2+qb,g*X3*X4*X2b+gb);




Oout[0]=Oout[0]+Ohelp;
/*****************************Left******************************/
if(BoundryLeft==1){

Ohelp.zeros(X0*X1*X2,X0*X1*X2);


for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Oascend[NumInterLeft](ib*X1+jb,i*X1+j);
Ohelp(ib*X2*X1+jb*X2+m,q)=sum;
}


Ohelp1.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int m=0;m<X2;m++)
sum=sum+Ohelp(i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+m,qb));
Ohelp1(qb,q)=sum;
}
Ohelp.eye(X2b*X3*X4,X2b*X3*X4);
Oout[0]=Oout[0]+kron(Ohelp1,Ohelp);
}
/*****************************Right******************************/
if(BoundryRight==1){


Ohelp.zeros(X3*X4*X2b,X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Oascend[NumInterRight](vb*X4+zb,v*X4+z);
Ohelp(x*X3*X4+vb*X4+zb,g)=sum;
}


Ohelp1.zeros(X2b*X3*X4,X2b*X3*X4);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
for(int x=0;x<X2b;x++)
sum=sum+Ohelp(x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](x*X3*X4+v*X4+z,gb));
Ohelp1(gb,g)=sum;
}
Ohelp.eye(X0*X1*X2,X0*X1*X2);

Oout[0]=Oout[0]+kron(Ohelp,Ohelp1);


}


}













