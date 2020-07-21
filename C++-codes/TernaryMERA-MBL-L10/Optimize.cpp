#include"Header.h"
//void Optimize(cx_mat*I0, cx_mat*I1, cx_mat*I2, cx_mat*I3, cx_mat*Itop, O*oH,O*oV,cx_mat*oF,double&Energy)
//{
//cx_mat * o1=new cx_mat[1];
//cx_mat * o2=new cx_mat[1];
//
//cx_mat * OM=new cx_mat[1];
//cx_mat * Or1=new cx_mat[Int];
//cx_mat * Od1=new cx_mat[Int];
//cx_mat * Ol1=new cx_mat[Int];
//cx_mat * Ou1=new cx_mat[Int];
//
//
//cx_mat * Or2=new cx_mat[Int];
//cx_mat * Od2=new cx_mat[Int];
//cx_mat * Ol2=new cx_mat[Int];
//cx_mat * Ou2=new cx_mat[Int];
//
//
//cx_mat * OHelp=new cx_mat[1];
//cx_mat * OHelp1=new cx_mat[1];
//cx_mat * OHelp2=new cx_mat[1];
//OM[0].zeros(Xi[1],Xi[1]);
//
//for(int q=0; q<Int; q++){
//Or1[q].zeros(Xi[1],Xi[1]);
//Ol1[q].zeros(Xi[1],Xi[1]);
//Ou1[q].zeros(Xi[1],Xi[1]);
//Od1[q].zeros(Xi[1],Xi[1]);
//Or2[q].zeros(Xi[1],Xi[1]);
//Ol2[q].zeros(Xi[1],Xi[1]);
//Ou2[q].zeros(Xi[1],Xi[1]);
//Od2[q].zeros(Xi[1],Xi[1]);
//}
//OHelp[0].zeros(Xi[1],Xi[1]);
//OHelp1[0].zeros(Xi[1],Xi[1]);
//OHelp2[0].zeros(Xi[1],Xi[1]);
//
//int iI;
//int iS;
//int jS;
//int jI;
///*****************/
//AscendingF(oF,I1,OHelp);
//iS=1;
//AscendingT1site( OHelp, Itop, OM,iS);
//
//AscendingF(oF,I2,OHelp);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//OM[0]=OM[0]+OHelp1[0];
//AscendingF(oF,I3,OHelp);
//iS=3;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//OM[0]=OM[0]+OHelp1[0];
///******************/
//Ascendingfull(oH,oV,I1,OHelp);
//iS=1;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//OM[0]=OM[0]+OHelp1[0];
//Ascendingfull(oH,oV,I2,OHelp);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//OM[0]=OM[0]+OHelp1[0];
//
//Ascendingfull(oH,oV,I3,OHelp);
//iS=3;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//OM[0]=OM[0]+OHelp1[0];
///*******************/
//for(int q=0; q<Int; q++){
//iS=2;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=0;
//iI=1;
//o2[0]=oV[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=1;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
//
///**/
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=1;
//iI=1;
//o2[0]=oV[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=1;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
//
//
///*******************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=2;
//iI=0;
//o2[0]=oV[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=1;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
//
///**/
//for(int q=0; q<Int; q++){
//iS=1;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=3;
//iI=0;
//o2[0]=oV[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=1;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
///********************************************/
//for(int q=0; q<Int; q++){
//iS=1;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=0;
//iI=1;
//o2[0]=oH[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=2;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
///****/
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=2;
//iI=1;
//o2[0]=oH[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=2;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
//
///*******************************************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=1;
//iI=0;
//o2[0]=oH[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=2;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
///****/
//for(int q=0; q<Int; q++){
//iS=2;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=3;
//iI=0;
//o2[0]=oH[0].M[q][iI];
//Ascending1site( o2, I3, OHelp1, iS);
//iS=2;
//jS=3;
//AscendingT2site(OHelp,OHelp1, Itop, OHelp2,iS,jS);
//OM[0]=OM[0]+OHelp2[0];}
//
//
///**************************************************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=1;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Or1[q]=OHelp1[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=2;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=1;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Or2[q]=OHelp1[0];
//}
///********************************************************/
//for(int q=0; q<Int; q++){
//iS=1;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=1;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Ol1[q]=OHelp1[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//Ascending1site( o1, I1, OHelp, iS);
//iS=1;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Ol2[q]=OHelp1[0];
//}
///********************************************************/
//for(int q=0; q<Int; q++){
//iS=2;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Ou1[q]=OHelp1[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Ou2[q]=OHelp1[0];
//}
///********************************************************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Od1[q]=OHelp1[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=1;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//Ascending1site( o1, I2, OHelp, iS);
//iS=2;
//AscendingT1site( OHelp, Itop, OHelp1,iS);
//Od2[q]=OHelp1[0];
//}
///********************************************************/
//
//int D=Xi[1]*Xi[0]*Xi[0]*Xi[0]*Xi[0];
//int Xl=D;
//int m=10;
//int W=2;
//
///////////////////////////////////////
//int num=0;
//double E1;
//cx_vec *Vechelp=new cx_vec[1];
//Vechelp[0]=zeros<cx_vec>(D);
//cx_vec * q;   //construct krylof's space
//cx_vec  a(m+1);              //elements of matrix of T=trans(Q)AQ
//cx_vec  b(m+1);              //elements of matrix of T=trans(Q)AQ
//cx_vec Vec;                  //final vectore of A*Vec=a*Vec
//cx_vec Vec1;                  //final vectore of A*Vec=a*Vec
//vec eigval1;                 //eigenval of T that is also eigenval of A
//cx_mat eigvec1;
//cx_mat V(m,m);
//cx_vec* r=new cx_vec[1];
//r[0]=randn<cx_vec>(D);
//cx_mat Q;
//
//for(int i=0; i<Xi[0];i++)
// for(int j=0; j<Xi[0];j++)
//   for(int m=0; m<Xi[0];m++)
//     for(int n=0; n<Xi[0];n++)
//       for(int p=0; p<Xi[1];p++)
//r[0](i*Xi[0]*Xi[0]*Xi[0]*Xi[1]+j*Xi[0]*Xi[0]*Xi[1]+m*Xi[0]*Xi[1]+n*Xi[1]+p)=I0[1](i*Xi[0]*Xi[0]*Xi[0]+j*Xi[0]*Xi[0]+m*Xi[0]+n,p);
//Vec=r[0]+0.000001*randn<cx_vec>(D);
//
//for(int p=0;p<=W;p++)
//{
//
//    if(p==0)
//q=new cx_vec[m+1];
//
//for(int j=0;j<m+1;j++)
//q[j].zeros(D);
//a.zeros(m+1);
//b.zeros(m+1);
//r[0]=Vec;
//b[0]=norm(r[0],2);
//for(int i=0;i<m;i++){
//if(norm(r[0],2)>1.0e-7)
//q[i+1]=r[0]/(b(i));
//else
//break;
//Vechelp[0]=q[i+1];
//Multi(oH,oV,oF,OM,Or1,Or2,Ol1,Ol2,Od1,Od2,Ou1,Ou2,Vechelp,D);
//
//r[0]=(Vechelp[0])-(b(i)*q[i]);
//
//
//a(i+1)=trace(trans(q[i+1])*r[0]);
//r[0]=r[0]-(a(i+1)*q[i+1]);
//
//
////orthogonalization
//{
//for(int n=1;n<=i;n++)//r is determined as to be orthonormal to all other vectors<=i
//r[0]=r[0]-(trace(trans(q[n])*r[0])*q[n]);
//}
//if(i<=m-2)
//b(i+1)=norm(r[0],2);
//}
//
//V.zeros(m,m);
//
//for(int i=1;i<m-1;i++)
//{
//V(i,i-1)=b(i);
// V(i,i)=a(i+1);
//  V(i,i+1)=b(i+1);
//}
//V(0,0)=a(1);
//V(0,1)=b(1);
//V(m-1,m-1)=a(m);
//V(m-1,m-2)=b(m-1);
//
//
//eig_sym(eigval1, eigvec1, V, "standard");
//
//Q.zeros(D,m);
//for(int n=0;n<m;n++)
//     Q.col(n)=q[n+1];
//       Vec=Q*eigvec1.col(0);
//
//
//
//if(p==W && num==0){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//num++;
//delete [] q;
//q=NULL;
//}else if(p==W)
//{
//    num++;
//    if(abs(eigval1(0))>0.000001){
//if(   (((abs(eigval1(0)-E1))/(abs(eigval1(0))))< 0.000000001)  )
//num++;
//else if(m<=50){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//delete [] q;
//q=NULL;
//}}else{
//if(  ((abs(eigval1(0)-E1))< 0.0000000001)  )
//num++;
//else if(m<=50){
//p=-1;
//m=m+10;
//E1=eigval1(0);
//delete [] q;
//q=NULL;
//}
//}
//}
//
//
//}
//
//
//cout<< setprecision(15)<<"  Energy="<<
//((eigval1(0))+Energy)/64.00<<" m="<<m<<endl;
//
//
//Vec=Vec/norm(Vec,2);
//for(int i=0; i<Xi[0];i++)
// for(int j=0; j<Xi[0];j++)
//   for(int m=0; m<Xi[0];m++)
//     for(int n=0; n<Xi[0];n++)
//       for(int p=0; p<Xi[1];p++)
//I0[0](i*Xi[0]*Xi[0]*Xi[0]+j*Xi[0]*Xi[0]+m*Xi[0]+n,p)=Vec(i*Xi[0]*Xi[0]*Xi[0]*Xi[1]+j*Xi[0]*Xi[0]*Xi[1]+m*Xi[0]*Xi[1]+n*Xi[1]+p);
//I0[1]=I0[0];
//
//
//
//delete [] q;
//delete [] r;
//delete [] Vechelp;
//delete [] o1;delete [] o2;
//delete [] OM;delete [] Or1;delete [] Or2;
//delete [] Od1;delete [] Od2;delete [] Ol1;delete [] Ol2;
//delete [] Ou1;delete [] Ou2;delete [] OHelp;
//delete [] OHelp1;delete [] OHelp2;
//}
//
//
//
//void Multi( O*oH,O*oV,cx_mat*oF,cx_mat*OM,cx_mat*Or1,cx_mat*Or2,cx_mat*Ol1,cx_mat*Ol2,
//cx_mat*Od1,cx_mat*Od2,cx_mat*Ou1,cx_mat*Ou2,cx_vec*VecHelp, int&D )
//{
//cx_mat * o1=new cx_mat[1];
//cx_mat * o2=new cx_mat[1];
//cx_mat *I=new cx_mat[1];
//cx_mat *I1=new cx_mat[1];
//cx_mat *Iresult=new cx_mat[1];
//cx_mat * Oact4=new cx_mat[1];
//I[0].zeros(Xi[0]*Xi[0]*Xi[0]*Xi[0],Xi[1]);
//I1[0].zeros(Xi[0]*Xi[0]*Xi[0]*Xi[0],Xi[1]);
//Iresult[0].zeros(Xi[0]*Xi[0]*Xi[0]*Xi[0],Xi[1]);
//for(int i=0; i<Xi[0];i++)
// for(int j=0; j<Xi[0];j++)
//   for(int m=0; m<Xi[0];m++)
//     for(int n=0; n<Xi[0];n++)
//       for(int p=0; p<Xi[1];p++)
//I[0](i*Xi[0]*Xi[0]*Xi[0]+j*Xi[0]*Xi[0]+m*Xi[0]+n,p)=VecHelp[0](i*Xi[0]*Xi[0]*Xi[0]*Xi[1]+j*Xi[0]*Xi[0]*Xi[1]+m*Xi[0]*Xi[1]+n*Xi[1]+p);
//I1[0]=I[0];
//
//
//RefreshF(oF,I1);
//Iresult[0]=I1[0]+Iresult[0];
///************************************/
//I1[0]=I[0];
//int iS=0;
//int iI=0;
//int jS=1;
//int jI=1;
//Oact4[0].eye(Xi[1],Xi[1]);
//int S4;
//S4=0;
//Refresh4site( oH,oV, Oact4, I1,S4);
//Iresult[0]=I1[0]+Iresult[0];
///****************************************************************/
//I1[0]=I[0];
//I1[0]=I1[0]*OM[0];
//Iresult[0]=I1[0]+Iresult[0];
///*****************************************************************/
//for(int q=0; q<Int; q++){
//iS=1;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//o2[0]=Or1[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oH[0].M[q][iI];
//o2[0]=Or2[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
///********************************************************************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//o2[0]=Ol1[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=2;
//iI=1;
//o1[0]=oH[0].M[q][iI];
//o2[0]=Ol2[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
//
///********************************************************************/
//for(int q=0; q<Int; q++){
//iS=2;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//o2[0]=Od1[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=3;
//iI=0;
//o1[0]=oV[0].M[q][iI];
//o2[0]=Od2[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
///********************************************************************/
//for(int q=0; q<Int; q++){
//iS=0;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//o2[0]=Ou1[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
//
//for(int q=0; q<Int; q++){
//iS=1;
//iI=1;
//o1[0]=oV[0].M[q][iI];
//o2[0]=Ou2[q];
// S4=1;
// I1[0]=I[0];
// Refresh1site( o1, o2, I1,  iS, S4);
//Iresult[0]=I1[0]+Iresult[0];
//}
///********************************************************************/
//for(int i=0; i<Xi[0];i++)
// for(int j=0; j<Xi[0];j++)
//   for(int m=0; m<Xi[0];m++)
//     for(int n=0; n<Xi[0];n++)
//       for(int p=0; p<Xi[1];p++)
//VecHelp[0](i*Xi[0]*Xi[0]*Xi[0]*Xi[1]+j*Xi[0]*Xi[0]*Xi[1]+m*Xi[0]*Xi[1]+n*Xi[1]+p)=Iresult[0](i*Xi[0]*Xi[0]*Xi[0]+j*Xi[0]*Xi[0]+m*Xi[0]+n,p);
//delete [] Iresult;
//delete [] I;
//delete [] I1;
//delete [] Oact4;
//delete [] o1;
//delete [] o2;
//}

