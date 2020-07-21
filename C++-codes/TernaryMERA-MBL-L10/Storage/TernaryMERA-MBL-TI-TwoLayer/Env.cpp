#include"Header.h"
void EnvIsoLeft(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,int BoundryRight,int BoundryLeft,int NumInterLeft,
int NumInterRight,vec Xlvec,cx_mat* Oout)
{

wall_clock timer;
double n_secs;

    complex<double> sum;
    cx_mat Ohelp;
    cx_mat Ohelp1;
    cx_mat Ohelp2;
int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);


/***************Interaction-Midle*********************************/
//timer.tic();

Ohelp1.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp1(g,xb*X2b+x)=sum;
}
Ohelp2=Rho[NumRhoup]*Ohelp1;


Ohelp=trans(Uni[NumUni])*Oascend[NumInter+1]*Uni[NumUni];

Ohelp1.zeros(X2b*X2b,X2*X2);
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b+x,mb*X2+m)=
Ohelp(mb*X2b+xb,m*X2b+x);

Ohelp=Ohelp2*(Ohelp1);




Ohelp1.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++){sum=0;
for(int mb=0;mb<X2;mb++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q))*Ohelp(q,mb*X2+m);
Ohelp1(q,i*X2*X1+j*X2+m)=sum;
}



Oout[0]=Ohelp1;
//cout<<trace(Ohelp1*Isometry[NumIso])<<endl;
//n_secs = timer.toc();
//cout << "took0 " << n_secs << " seconds" << endl;

/***************************Interaction-Left*********************/
//timer.tic();

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp2=Rho[NumRhoup]*Ohelp;



Ohelp.zeros(X2b,X2*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;
//timer.tic();

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
//n_secs = timer.toc();
//cout << "took2 " << n_secs << " seconds" << endl;
//timer.tic();

Ohelp1.zeros(X2b*X2b,X2*X2*X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b+x,jb*X1*X2*X2+j*X2*X2+mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j);


Ohelp=Ohelp2*Ohelp1;
//n_secs = timer.toc();
//cout << "took3 " << n_secs << " seconds" << endl;
//timer.tic();
Ohelp1.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++){sum=0;
for(int mb=0;mb<X2;mb++)
for(int jb=0;jb<X1;jb++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,q))*Ohelp(q,jb*X1*X2*X2+j*X2*X2+mb*X2+m);
Ohelp1(q,i*X2*X1+j*X2+m)=sum;
}

Oout[0]=Ohelp1+Oout[0];
//n_secs = timer.toc();
//cout << "took4 " << n_secs << " seconds" << endl;
//timer.tic();

//cout<<trace(Oout[0]*Isometry[NumIso])<<endl;

//n_secs = timer.toc();
//cout << "took5 " << n_secs << " seconds" << endl;
/****************************Boundry-Right*************************************/
if(BoundryRight==1){

Ohelp2.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
sum=sum+conj(Isometry[NumIso+1](x*X3*X4+v*X4+z,g))*Oascend[NumInterRight](v*X4+z,vb*X4+zb);
Ohelp2(x*X3*X4+vb*X4+zb,g)=sum;
}


Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Ohelp2(xb*X4*X3+v*X4+z,g);
Ohelp(g,xb*X2b+x)=sum;
}
Ohelp2=Rho[NumRhoup]*Ohelp;



Ohelp.zeros(X2b,X2*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(c,s*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;
//timer.tic();

Ohelp1.zeros(X2b*X2b*X2*X2,X2*X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,sb*X2+s)=Ohelp(sb*X2b*X2+mb*X2b+xb,s*X2b*X2+m*X2b+x);


cx_mat OascendH;
OascendH.eye(X2*X1,X2*X1);

Ohelp.eye(X2*X2,X1*X1);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(sb*X2+s,jb*X1+j)=OascendH(jb*X2+sb,j*X2+s);


Ohelp=Ohelp1*Ohelp;
//n_secs = timer.toc();
//cout << "took2 " << n_secs << " seconds" << endl;
//timer.tic();

Ohelp1.zeros(X2b*X2b,X2*X2*X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(xb*X2b+x,jb*X1*X2*X2+j*X2*X2+mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j);


Ohelp=Ohelp2*Ohelp1;
//n_secs = timer.toc();
//cout << "took3 " << n_secs << " seconds" << endl;
//timer.tic();
Ohelp1.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++){sum=0;
for(int mb=0;mb<X2;mb++)
for(int jb=0;jb<X1;jb++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,q))*Ohelp(q,jb*X1*X2*X2+j*X2*X2+mb*X2+m);
Ohelp1(q,i*X2*X1+j*X2+m)=sum;
}

Oout[0]=Ohelp1+Oout[0];



}
/***************************Interaction-Right*********************/
//timer.tic();

Ohelp.zeros(X2,X2b*X2*X2b);
for(int s=0;s<X2;s++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(s,c*X2b*X2+m*X2b+x)=Uni[NumUni](s*X2b+c,m*X2b+x);

Ohelp=trans(Ohelp)*Ohelp;
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;

//timer.tic();
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
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;

//timer.tic();

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

Ohelp2=Ohelp*Ohelp2;
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;
//timer.tic();

Ohelp=Rho[NumRhoup]*Ohelp2;



Ohelp1.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++){sum=0;
for(int mb=0;mb<X2;mb++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q))*Ohelp(q,mb*X2+m);
Ohelp1(q,i*X2*X1+j*X2+m)=sum;
}

Oout[0]=Ohelp1+Oout[0];

//cout<<trace(Oout[0]*Isometry[NumIso])<<endl;
//n_secs = timer.toc();
//cout << "took1 " << n_secs << " seconds" << endl;

/******************************************Boundry-Left***************************************************/
//timer.tic();

if(BoundryLeft==1){


Ohelp.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+j*X2+m,q))*Oascend[NumInterLeft](i*X1+j,ib*X1+jb);
Ohelp(ib*X2*X1+jb*X2+m,q)=sum;
}

Ohelp1.zeros(X0*X1*X2,1);
for(int q=0;q<X0*X1*X2;q++){sum=0;
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Rho[NumRhoup](q,g);
Ohelp1(q,0)=sum;
}

Ohelp2.zeros(X2*X1*X0,X2*X1*X0);
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
Ohelp2(q,i*X2*X1+j*X2+m)=Ohelp1(q,0)*Ohelp(i*X2*X1+j*X2+m,q);


Oout[0]=Ohelp2+Oout[0];


}
//n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;

}

void EnvIsoRight(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,int BoundryRight,int BoundryLeft,int NumInterLeft,
int NumInterRight,vec Xlvec,cx_mat* Oout)
{
//wall_clock timer;
//double n_secs;

complex<double> sum;
cx_mat Ohelp;
cx_mat Ohelp1;
cx_mat Ohelp2;
int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);

/*********************B-Midle:Interaction-Right**********************/
//timer.tic();

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

Ohelp1=strans(Rho[NumRhoup])*Ohelp1;


Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++){sum=0;
for(int xb=0;xb<X2b;xb++)
sum=sum+conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g))*Ohelp1(g,xb*X2b+x);
Ohelp(g,x*X4*X3+v*X4+z)=sum;
}


Oout[0]=Ohelp;
//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
//n_secs = timer.toc();

//cout << "took " << n_secs << " seconds" << endl;

/*******************************Interaction-Left*************************/
//timer.tic();

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


Ohelp1=strans(Rho[NumRhoup])*Ohelp2;


Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++){sum=0;
for(int xb=0;xb<X2b;xb++)
sum=sum+conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g))*Ohelp1(g,xb*X2b+x);
Ohelp(g,x*X4*X3+v*X4+z)=sum;
}


Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
//n_secs = timer.toc();

//cout << "took " << n_secs << " seconds" << endl;

/************************Interaction-Right**********************************************/
//timer.tic();

Ohelp.zeros(X2*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp2=Ohelp*Rho[NumRhoup];


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

Ohelp1.zeros(X2b*X2b*X3*X3,X2*X2);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v);


Ohelp1=Ohelp1*Ohelp2;


Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++){sum=0;
for(int xb=0;xb<X2b;xb++)
for(int vb=0;vb<X3;vb++)
sum=sum+conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,g))*Ohelp1(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,g);
Ohelp(g,x*X4*X3+v*X4+z)=sum;
}


Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
//n_secs = timer.toc();

//cout << "took " << n_secs << " seconds" << endl;

/***********************BoundryRight***********************************************/
if(BoundryLeft==1){


Ohelp2.zeros(X2*X1*X0,X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+conj(Isometry[NumIso](i*X2*X1+j*X2+m,q))*Oascend[NumInterLeft](i*X1+j,ib*X1+jb);
Ohelp2(ib*X2*X1+jb*X2+m,q)=sum;
}



Ohelp.zeros(X2*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Ohelp2(i*X2*X1+j*X2+mb,q);
Ohelp(mb*X2+m,q)=sum;
}

Ohelp2=Ohelp*Rho[NumRhoup];


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



cx_mat OascendH;
OascendH.eye(X3*X2b,X3*X2b);
Ohelp.zeros(X2b*X2b,X3*X3);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(cb*X2b+c,vb*X3+v)=OascendH(cb*X3+vb,c*X3+v);


Ohelp=Ohelp1*Ohelp;

Ohelp1.zeros(X2b*X2b*X3*X3,X2*X2);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,mb*X2+m)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,vb*X3+v);


Ohelp1=Ohelp1*Ohelp2;


Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++){sum=0;
for(int xb=0;xb<X2b;xb++)
for(int vb=0;vb<X3;vb++)
sum=sum+conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,g))*Ohelp1(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,g);
Ohelp(g,x*X4*X3+v*X4+z)=sum;
}


Oout[0]=Ohelp+Oout[0];




}



//timer.tic();

if(BoundryRight==1){


Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
sum=sum+conj(Isometry[NumIso+1](x*X3*X4+v*X4+z,g))*Oascend[NumInterRight](v*X4+z,vb*X4+zb);
Ohelp(x*X3*X4+vb*X4+zb,g)=sum;
}


Ohelp1.zeros(X3*X4*X2b,1);
for(int g=0;g<X3*X4*X2b;g++){sum=0;
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Rho[NumRhoup](q,g);
Ohelp1(g,0)=sum;
}

Ohelp2.zeros(X3*X4*X2b,X3*X4*X2b);
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
for(int x=0;x<X2b;x++)
for(int g=0;g<X3*X4*X2b;g++)
Ohelp2(g,x*X3*X4+v*X4+z)=Ohelp1(g,0)*Ohelp(x*X3*X4+v*X4+z,g);

Oout[0]=Ohelp2+Oout[0];
}
//n_secs = timer.toc();

//cout << "took " << n_secs << " seconds" << endl;


}

void EnvUni(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,int BoundryRight,int BoundryLeft,int NumInterLeft,int NumInterRight
,vec Xlvec,cx_mat* Oout)
{

complex<double> sum;
cx_mat Ohelp;
cx_mat Ohelp1;
cx_mat Ohelp2;
int X0=Xlvec(0);
int X1=Xlvec(1);
int X2=Xlvec(2);
int X2b=Xlvec(3);
int X3=Xlvec(4);
int X4=Xlvec(5);
/************************Midle************************/
Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Ohelp1=Rho[NumRhoup]*Ohelp;


Ohelp.zeros(X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(q,mb*X2+m)=sum;
}



Ohelp1=strans(Ohelp)*Ohelp1;

Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp1(mb*X2+m,xb*X2b+x);


Ohelp=Ohelp*trans(Uni[NumUni])*Oascend[NumInter+1];

Oout[0]=Ohelp;
//cout<<trace(Oout[0]*Uni[NumUni])<<endl;
/*******************************Boundry-Left*******************************************/

if(BoundryLeft==1){




Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Ohelp1=Rho[NumRhoup]*Ohelp;

Ohelp2.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int ib=0;ib<X0;ib++)
for(int jb=0;jb<X1;jb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*Oascend[NumInterLeft](ib*X1+jb,i*X1+j);
Ohelp2(ib*X2*X1+jb*X2+m,q)=sum;
}


Ohelp.zeros(X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Ohelp2(i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(q,mb*X2+m)=sum;
}



Ohelp1=strans(Ohelp)*Ohelp1;

Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp1(mb*X2+m,xb*X2b+x);


Ohelp=Ohelp*trans(Uni[NumUni]);

Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Uni[NumUni])<<endl;
}
/********************************Right-Boundry*********************************************************/
if(BoundryRight==1){

Ohelp2.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int vb=0;vb<X3;vb++)
for(int zb=0;zb<X4;zb++){sum=0;
for(int z=0;z<X4;z++)
for(int v=0;v<X3;v++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*Oascend[NumInterRight](vb*X4+zb,v*X4+z);
Ohelp2(x*X3*X4+vb*X4+zb,g)=sum;
}



Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Ohelp2(x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Ohelp1=Rho[NumRhoup]*Ohelp;



Ohelp.zeros(X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(q,mb*X2+m)=sum;
}



Ohelp1=strans(Ohelp)*Ohelp1;

Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp1(mb*X2+m,xb*X2b+x);


Ohelp=Ohelp*trans(Uni[NumUni]);

Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Uni[NumUni])<<endl;
}




/*************************Left***************************************/

Ohelp.zeros(X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,g));
Ohelp(g,xb*X2b+x)=sum;
}

Ohelp1=Rho[NumRhoup]*Ohelp;


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

Ohelp1.zeros(X2*X1*X1*X2b,X2*X2b);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(m*X2b*X1*X1+x*X1*X1+jb*X1+j,mb*X2b+xb)=
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x);

Ohelp2=Ohelp1*trans(Uni[NumUni]);


Ohelp1.zeros(X2*X2b*X2b,X2*X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int sb=0;sb<X2;sb++)
Ohelp1(m*X2b*X2b+x*X2b+cb,jb*X2*X1+j*X2+sb)=
Ohelp2(m*X2b*X1*X1+x*X1*X1+jb*X1+j,sb*X2b+cb);


Ohelp.zeros(X2*X1*X1,X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(jb*X2*X1+j*X2+sb,s)=Oascend[NumInter](jb*X2+sb,j*X2+s);


Ohelp1=Ohelp1*Ohelp;


Ohelp.zeros(X2b*X2,X2b*X2);
for(int s=0;s<X2;s++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,s*X2b+cb)=Ohelp1(m*X2b*X2b+x*X2b+cb,s);


Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Uni[NumUni])<<endl;
/******************************************************************************/

Ohelp.zeros(X2*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,q));
Ohelp(mb*X2+m,q)=sum;
}

Ohelp1=Ohelp*Rho[NumRhoup];

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

Ohelp1.zeros(X2*X3*X3*X2b,X2*X2b);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(m*X2b*X3*X3+x*X3*X3+vb*X3+v,mb*X2b+xb)=
Ohelp2(mb*X2+m,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x);


Ohelp2=Ohelp1*trans(Uni[NumUni]);


Ohelp1.zeros(X2*X2*X2b,X2b*X3*X3);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int sb=0;sb<X2;sb++)
Ohelp1(m*X2b*X2+x*X2+sb,vb*X2b*X3+v*X2b+cb)=
Ohelp2(m*X2b*X3*X3+x*X3*X3+vb*X3+v,sb*X2b+cb);


Ohelp.zeros(X3*X2b*X3,X2b);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(vb*X3*X2b+v*X2b+cb,c)=Oascend[NumInter+2](cb*X3+vb,c*X3+v);

Ohelp1=Ohelp1*Ohelp;


Ohelp.zeros(X2b*X2,X2b*X2);
for(int sb=0;sb<X2;sb++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,sb*X2b+c)=Ohelp1(m*X2b*X2+x*X2+sb,c);


Oout[0]=Ohelp+Oout[0];
}


void EnvIsoLeftMERA(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,int BoundryLeft,
int NumInterLeft,vec Xlvec,cx_mat* Oout){
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

/*****************************Midle**************************************************/
Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(gb*X3*X4*X2b+g,xb*X2b+x)=sum;
}
Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp=Ohelp1*Ohelp;

Ohelp2=trans(Uni[NumUni])*Oascend[NumInter+1]*Uni[NumUni];

Ohelp1.zeros(X2b*X2b,X2*X2);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp1(xb*X2b+x,mb*X2+m)=
Ohelp2(mb*X2b+xb,m*X2b+x);

Ohelp=Ohelp*Ohelp1;


Ohelp1.zeros(X0*X1*X2*X2,X0*X1*X2*X2);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp1(m*X0*X1*X2+q,qb*X2+mb)=
Ohelp(qb*X0*X1*X2+q,mb*X2+m);


Ohelp.zeros(X0*X1*X2*X2,X0*X1);
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int mb=0;mb<X2;mb++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp(qb*X2+mb,i*X1+j)=
conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));

Ohelp1=Ohelp1*Ohelp;

Ohelp.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++)
Ohelp(q,i*X1*X2+j*X2+m)=Ohelp1(m*X0*X1*X2+q,i*X1+j);

Oout[0]=Ohelp;
//cout<<endl<<endl;
//cout<<trace(Ohelp*Isometry[NumIso])<<endl;

/*********************Left**********************************************/

Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(gb*X3*X4*X2b+g,xb*X2b+x)=sum;
}

Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp2=Ohelp1*Ohelp;

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
Ohelp1.zeros(X2*X2*X1*X1,X2b*X2b);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp1(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x)=Ohelp(xb*X2b*X2*X2+x*X2*X2+mb*X2+m,jb*X1+j);

Ohelp=Ohelp1*strans(Ohelp2);

Ohelp1.zeros(X0*X1*X2*X2*X1,X0*X1*X2*X2*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp1(j*X0*X1*X2*X2+m*X0*X1*X2+q,jb*X0*X1*X2*X2+mb*X0*X1*X2+qb)=Ohelp(jb*X1*X2*X2+j*X2*X2+mb*X2+m,qb*X0*X1*X2+q);

Ohelp.zeros(X0*X1*X2*X2*X1,X0);
for(int i=0;i<X0;i++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++)
for(int qb=0;qb<X0*X1*X2;qb++)
Ohelp(jb*X0*X1*X2*X2+mb*X0*X1*X2+qb,i)=
conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,qb));

Ohelp1=Ohelp1*Ohelp;

Ohelp.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++)
Ohelp(q,i*X1*X2+j*X2+m)=Ohelp1(j*X0*X1*X2*X2+m*X0*X1*X2+q,i);

Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Isometry[NumIso])<<endl;
/******************************Right************************************************************/

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
Ohelp(gb*X3*X4*X2b+g,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x)=sum;
}

Ohelp2=Ohelp*Ohelp2;



Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);


Ohelp1=Ohelp1*Ohelp2;


Ohelp.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int i=0;i<X0;i++){sum=0;
for(int qb=0;qb<X0*X1*X2;qb++)
for(int mb=0;mb<X2;mb++)
sum=sum+Ohelp1(qb*X0*X1*X2+q,mb*X2+m)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(q,i*X1*X2+j*X2+m)=sum;
}


Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Isometry[NumIso])<<endl;
/***********************************************************/
if(BoundryLeft==1){

cx_mat Identitiy;
Identitiy.eye(X2,X2);

Ohelp=kron(Oascend[NumInterLeft],Identitiy);
Ohelp=trans(Isometry[NumIso])*Ohelp;

Ohelp1.zeros(X0*X1*X2,X0*X1*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++){sum=0;
for(int g=0;g<X3*X4*X2b;g++)
sum=sum+Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+g);
Ohelp1(qb,q)=sum;
}

Ohelp2.zeros(X2*X1*X0,X2*X1*X0);
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
for(int q=0;q<X0*X1*X2;q++)
for(int m=0;m<X2;m++){sum=0;
for(int qb=0;qb<X0*X1*X2;qb++)
sum=sum+Ohelp1(qb,q)*Ohelp(qb,i*X2*X1+j*X2+m);
Ohelp2(q,i*X2*X1+j*X2+m)=sum;}


Oout[0]=Ohelp2+Oout[0];

}


}



void EnvIsoRightMERA(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,int BoundryRight,
int NumInterRight,vec Xlvec,cx_mat* Oout){
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
/*******************************************Midle***********************************************************/
Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(qb*X0*X1*X2+q,mb*X2+m)=sum;
}
Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp=strans(Ohelp1)*Ohelp;


Ohelp2=trans(Uni[NumUni])*Oascend[NumInter+1]*Uni[NumUni];

Ohelp1.zeros(X2*X2,X2b*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp1(mb*X2+m,xb*X2b+x)=
Ohelp2(mb*X2b+xb,m*X2b+x);

Ohelp1=Ohelp*Ohelp1;

Ohelp.zeros(X3*X4*X2b*X2b,X3*X4*X2b*X2b);
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp(x*X3*X4*X2b+g,gb*X2b+xb)=
Ohelp1(gb*X3*X4*X2b+g,xb*X2b+x);


Ohelp2.zeros(X3*X4*X2b*X2b,X3*X4);
for(int xb=0;xb<X2b;xb++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp2(gb*X2b+xb,v*X4+z)=
conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));

Ohelp1=Ohelp*Ohelp2;

Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int g=0;g<X3*X4*X2b;g++)
Ohelp(g,x*X4*X3+v*X4+z)=Ohelp1(x*X3*X4*X2b+g,v*X4+z);

Oout[0]=Ohelp;

//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
/********************************************Left**********************************/

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
Ohelp(qb*X2*X1*X0+q,jb*X1*X2*X2+j*X2*X2+mb*X2+m)=sum;
}
Ohelp2=Ohelp*Ohelp2;



Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);


Ohelp1=strans(Ohelp1)*Ohelp2;


Ohelp.zeros(X3*X4*X2b*X2b,X3*X4*X2b*X2b);
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp(x*X3*X4*X2b+g,gb*X2b+xb)=
Ohelp1(gb*X3*X4*X2b+g,xb*X2b+x);


Ohelp2.zeros(X3*X4*X2b*X2b,X3*X4);
for(int xb=0;xb<X2b;xb++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp2(gb*X2b+xb,v*X4+z)=
conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));

Ohelp1=Ohelp*Ohelp2;

Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int g=0;g<X3*X4*X2b;g++)
Ohelp(g,x*X4*X3+v*X4+z)=Ohelp1(x*X3*X4*X2b+g,v*X4+z);

Oout[0]=Ohelp+Oout[0];

//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
/**********************************************Right*********************************/


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

Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);

for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(qb*X0*X1*X2+q,mb*X2+m)=sum;
}

Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp=strans(Ohelp1)*Ohelp;

Ohelp1=Ohelp2*strans(Ohelp);


Ohelp2.zeros(X3*X4*X2b*X3*X2b,X3*X4*X2b*X3*X2b);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
Ohelp2(g*X2b*X3+x*X3+v,gb*X2b*X3+xb*X3+vb)=Ohelp1(vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x,gb*X3*X4*X2b+g);



Ohelp1.zeros(X3*X4*X2b*X3*X2b,X4);
for(int xb=0;xb<X2b;xb++)
for(int vb=0;vb<X3;vb++)
for(int z=0;z<X4;z++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(gb*X2b*X3+xb*X3+vb,z)=
conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,gb));

Ohelp1=Ohelp2*Ohelp1;

Ohelp.zeros(X3*X4*X2b,X3*X4*X2b);
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int g=0;g<X3*X4*X2b;g++)
Ohelp(g,x*X4*X3+v*X4+z)=Ohelp1(g*X2b*X3+x*X3+v,z);
Oout[0]=Ohelp+Oout[0];

//cout<<trace(Oout[0]*Isometry[NumIso+1])<<endl;
/***************************************************************/
if(BoundryRight==1){

cx_mat Identitiy;
Identitiy.eye(X2b,X2b);
Ohelp=kron(Identitiy,Oascend[NumInterRight]);
Ohelp=trans(Isometry[NumIso+1])*Ohelp;

Ohelp1.zeros(X3*X4*X2b,X3*X4*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++){sum=0;
for(int q=0;q<X0*X1*X2;q++)
sum=sum+Rho[NumRhoup](q*X3*X4*X2b+g,q*X3*X4*X2b+gb);
Ohelp1(gb,g)=sum;
}
Ohelp2.zeros(X3*X4*X2b,X3*X4*X2b);

for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
for(int g=0;g<X3*X4*X2b;g++)
for(int x=0;x<X2b;x++){sum=0;
for(int gb=0;gb<X3*X4*X2b;gb++)
sum=sum+Ohelp1(gb,g)*Ohelp(gb,x*X4*X3+v*X4+z);
Ohelp2(g,x*X4*X3+v*X4+z)=sum;}

Oout[0]=Ohelp2+Oout[0];

}
}

void EnvUniMERA(cx_mat * Oascend,cx_mat *Rho, cx_mat *Uni,cx_mat * Isometry,
int NumInter,int NumRhoup, int NumIso,int NumUni,vec Xlvec,cx_mat* Oout){

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
/*************************Midle**********************************/
Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(gb*X3*X4*X2b+g,xb*X2b+x)=sum;
}
Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp1=Ohelp1*Ohelp;


Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(qb*X0*X1*X2+q,mb*X2+m)=sum;
}



Ohelp2=strans(Ohelp)*Ohelp1;


Ohelp.zeros(X2*X2b,X2*X2b);
for(int mb=0;mb<X2;mb++)
for(int m=0;m<X2;m++)
for(int xb=0;xb<X2b;xb++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,mb*X2b+xb)=Ohelp2(mb*X2+m,xb*X2b+x);


Ohelp=Ohelp*trans(Uni[NumUni])*Oascend[NumInter+1];

Oout[0]=Ohelp;
/******************************Left***************************/

Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int v=0;v<X3;v++)
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+v*X4+z,gb));
Ohelp(gb*X3*X4*X2b+g,xb*X2b+x)=sum;
}

Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp2=Ohelp1*Ohelp;

Ohelp.zeros(X2*X2*X1*X1,X2*X1*X0*X2*X1*X0);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+jb*X2+mb,qb));
Ohelp(jb*X1*X2*X2+j*X2*X2+mb*X2+m,qb*X2*X1*X0+q)=sum;
}

Ohelp2=Ohelp*Ohelp2;

Ohelp1.zeros(X2*X1*X1*X2b,X2*X2b);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(m*X2b*X1*X1+x*X1*X1+jb*X1+j,mb*X2b+xb)=
Ohelp2(jb*X1*X2*X2+j*X2*X2+mb*X2+m,xb*X2b+x);

Ohelp2=Ohelp1*trans(Uni[NumUni]);


Ohelp1.zeros(X2*X2b*X2b,X2*X1*X1);
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
for(int x=0;x<X2b;x++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int sb=0;sb<X2;sb++)
Ohelp1(m*X2b*X2b+x*X2b+cb,jb*X2*X1+j*X2+sb)=
Ohelp2(m*X2b*X1*X1+x*X1*X1+jb*X1+j,sb*X2b+cb);


Ohelp.zeros(X2*X1*X1,X2);
for(int s=0;s<X2;s++)
for(int sb=0;sb<X2;sb++)
for(int j=0;j<X1;j++)
for(int jb=0;jb<X1;jb++)
Ohelp(jb*X2*X1+j*X2+sb,s)=Oascend[NumInter](jb*X2+sb,j*X2+s);


Ohelp1=Ohelp1*Ohelp;

Ohelp.zeros(X2b*X2,X2b*X2);
for(int s=0;s<X2;s++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,s*X2b+cb)=Ohelp1(m*X2b*X2b+x*X2b+cb,s);


Oout[0]=Ohelp+Oout[0];
//cout<<trace(Oout[0]*Uni[NumUni])<<endl;
/*******************************************************************/
Ohelp.zeros(X0*X1*X2*X0*X1*X2,X2*X2);

for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++){sum=0;
for(int i=0;i<X0;i++)
for(int j=0;j<X1;j++)
sum=sum+Isometry[NumIso](i*X2*X1+j*X2+m,q)*conj(Isometry[NumIso](i*X2*X1+j*X2+mb,qb));
Ohelp(qb*X0*X1*X2+q,mb*X2+m)=sum;
}

Ohelp1.zeros(X0*X1*X2*X0*X1*X2,X3*X4*X2b*X3*X4*X2b);
for(int q=0;q<X0*X1*X2;q++)
for(int qb=0;qb<X0*X1*X2;qb++)
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
Ohelp1(qb*X0*X1*X2+q,gb*X3*X4*X2b+g)=Rho[NumRhoup](q*X3*X4*X2b+g,qb*X3*X4*X2b+gb);

Ohelp1=strans(Ohelp1)*Ohelp;


Ohelp.zeros(X3*X4*X2b*X3*X4*X2b,X2b*X2b*X3*X3);
for(int g=0;g<X3*X4*X2b;g++)
for(int gb=0;gb<X3*X4*X2b;gb++)
for(int x=0;x<X2b;x++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int xb=0;xb<X2b;xb++){sum=0;
for(int z=0;z<X4;z++)
sum=sum+Isometry[NumIso+1](x*X3*X4+v*X4+z,g)*conj(Isometry[NumIso+1](xb*X4*X3+vb*X4+z,gb));
Ohelp(gb*X3*X4*X2b+g,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x)=sum;
}

Ohelp2=strans(Ohelp1)*Ohelp;

Ohelp1.zeros(X2*X3*X3*X2b,X2*X2b);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int xb=0;xb<X2b;xb++)
for(int m=0;m<X2;m++)
for(int mb=0;mb<X2;mb++)
Ohelp1(m*X2b*X3*X3+x*X3*X3+vb*X3+v,mb*X2b+xb)=
Ohelp2(mb*X2+m,vb*X3*X2b*X2b+v*X2b*X2b+xb*X2b+x);


Ohelp2=Ohelp1*trans(Uni[NumUni]);

Ohelp1.zeros(X2*X2*X2b,X2b*X3*X3);
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
for(int x=0;x<X2b;x++)
for(int cb=0;cb<X2b;cb++)
for(int m=0;m<X2;m++)
for(int sb=0;sb<X2;sb++)
Ohelp1(m*X2b*X2+x*X2+sb,vb*X2b*X3+v*X2b+cb)=
Ohelp2(m*X2b*X3*X3+x*X3*X3+vb*X3+v,sb*X2b+cb);

Ohelp.zeros(X3*X2b*X3,X2b);
for(int c=0;c<X2b;c++)
for(int cb=0;cb<X2b;cb++)
for(int v=0;v<X3;v++)
for(int vb=0;vb<X3;vb++)
Ohelp(vb*X3*X2b+v*X2b+cb,c)=Oascend[NumInter+2](cb*X3+vb,c*X3+v);


Ohelp1=Ohelp1*Ohelp;

Ohelp.zeros(X2b*X2,X2b*X2);
for(int sb=0;sb<X2;sb++)
for(int c=0;c<X2b;c++)
for(int m=0;m<X2;m++)
for(int x=0;x<X2b;x++)
Ohelp(m*X2b+x,sb*X2b+c)=Ohelp1(m*X2b*X2+x*X2+sb,c);

Oout[0]=Ohelp+Oout[0];
/**************************************/
}


