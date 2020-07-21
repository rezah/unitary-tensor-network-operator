//#include"Header.h"
//void RhoLanczos(cx_mat*Oascend,cx_mat*Uni,
//cx_mat*Rho,
//cx_mat*Isometry){
//
//int NumLay;
//
//int D=Xi[Lay]*Xi[Lay]*Xi[Lay]*Xi[Lay];
//int Xlh=Xi[Lay];
//int m=12;
//int W=2;
//int num=0;
//complex<double> max1;
//complex<double> E1;
///////////////////////////////////////
//cx_mat h;
//cx_vec eigval1;
//cx_mat eigvec1;
//cx_vec  r;                   //very well-known vector of Lanczos
//cx_mat Q;            //trans(Q)AQ=h
//cx_vec * q;  //Lanczos's Vectors
//cx_vec VecF;          //final Vector
//cx_vec Vechelp;
//cx_vec*Pointer=&Vechelp;
//uword  index;
//
//VecF.zeros(D);
//
// for(int i=0; i<Xlh;i++)
//     for(int j=0; j<Xlh;j++)
//     for(int m=0; m<Xlh;m++)
//     for(int n=0; n<Xlh;n++)
//VecF(i*Xlh*Xlh*Xlh+j*Xlh*Xlh+m*Xlh+n)=
//Rho[Lay](i*Xlh+j,m*Xlh+n);
//
//VecF=VecF+0.000001*randn<cx_vec>(D);
//for(int p=0;p<=W;p++){
//
//if(p==0)
//q=new cx_vec[m];
//
//for(int j=0;j<m;j++)
//q[j].zeros(D);
//r=VecF;
//q[0]=r/norm(r,2);
//h.zeros(m,m);
//Q.zeros(D,m);
//
//for(int j=0;j<m;j++){
//
//Vechelp=q[j];
//
//
//for(int i1=0; i1<Xlh;i1++)
//     for(int j=0; j<Xlh;j++)
//     for(int m=0; m<Xlh;m++)
//     for(int n=0; n<Xlh;n++)
//Rho[Lay](i1*Xlh+j,m*Xlh+n)=
//Vechelp(i1*Xlh*Xlh*Xlh+j*Xlh*Xlh+m*Xlh+n);
//
//
//NumLay=Lay;
//Dscen(Oascend, Uni, Rho,Isometry,NumLay);
//
//
//for(int i1=0; i1<Xlh;i1++)
//     for(int j=0; j<Xlh;j++)
//     for(int m=0; m<Xlh;m++)
//     for(int n=0; n<Xlh;n++)
//Vechelp(i1*Xlh*Xlh*Xlh+j*Xlh*Xlh+m*Xlh+n)=
//Rho[Lay-1](i1*Xlh+j,m*Xlh+n);
//
//
//r=Vechelp;
//
//for(int i=0;i<=j;i++){
//h(i,j)=trace(trans(q[i])*r);
//r=r-(h(i,j)*q[i]);
//}
//
//if(j<m-1){
//h(j+1,j)=norm(r,2);
//if(norm(r,2)>1.0e-8)
//q[j+1]=r/norm(r,2);
//else
//break;
//}
//
//}
//eig_gen(eigval1, eigvec1, h, 'r');
//max1=eigval1.max(index);
//
//
//for(int n=0;n<m;n++)
//    Q.col(n)=q[n];
//
//VecF=Q*eigvec1.col(index);
//if(p==W && num==0){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//num++;
//delete [] q;
//}else if(p==W)
//{
//    num++;
//    if(abs(eigval1(0))>0.000001){
//if(   (((abs(eigval1(0)-E1))/(abs(eigval1(0))))< 0.000000001)  )
//num++;
//else if(m<=60){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//delete [] q;
//}}else{
//if(  ((abs(eigval1(0)-E1))< 0.0000000001)  )
//num++;
//else if(m<=60){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//delete [] q;
//}
//}
//}
//}
//delete [] q;
//if(norm(VecF,2)>1.0e-20)
//VecF=VecF/norm(VecF,2);
//
//for(int i=0; i<Xlh;i++)
//     for(int j=0; j<Xlh;j++)
//     for(int m=0; m<Xlh;m++)
//     for(int n=0; n<Xlh;n++)
//Rho[Lay](i*Xlh+j,m*Xlh+n)=
//VecF(i*Xlh*Xlh*Xlh+j*Xlh*Xlh+m*Xlh+n);
//
//
//Rho[Lay]=(Rho[Lay]+trans(Rho[Lay]))/2.00;
//Rho[Lay]=Rho[Lay]*(1.00/abs(trace(Rho[Lay])));
//
//
//NumLay=Lay;
//Dscen(Oascend, Uni, Rho,Isometry,NumLay);
//double Sum=0;
//for(int i=0; i<Xlh;i++)
//     for(int j=0; j<Xlh;j++)
//     for(int m=0; m<Xlh;m++)
//     for(int n=0; n<Xlh;n++)
//Sum=Sum+abs(Rho[Lay](i*Xlh+j,m*Xlh+n)-Rho[Lay-1](i*Xlh+j,m*Xlh+n));
//if(Sum>0.000001){
//Sum=1-Sum;
//cout<<endl<<"Lanczos="<<m<<" EigEnvlue="<<norm(max1)<<" Similarity="<<Sum<<endl;
//}
//}
