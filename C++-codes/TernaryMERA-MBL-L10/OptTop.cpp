//#include"Header.h"
//void OptTop(cx_mat*Oascend,cx_mat*Isometry,cx_mat*Uni,cx_mat*Rho,cx_mat*EnvU,cx_mat*Env){
//int Ihelp;
//cx_mat Ohelp;
//
//cx_mat * UniHelp=new cx_mat[1];
//cx_mat * IsoHelp=new cx_mat[1];
//
//
//
//
//vec * Checker=new vec[1];
//Checker[0].zeros(4);
//Checker[0](1)=220.094;
//
//for(int q=0;q<=3;q++){
//Ihelp=Lay-2;
//Ascen(Oascend, Uni, Rho,Isometry,Ihelp);
//Ohelp=Oascend[Lay-1];
//
//for(int i1=0;i1<q;i1++){
//Ascen(Oascend, Uni, Rho,Isometry,Lay-1);
//Ohelp=Ohelp+Oascend[Lay];
//Oascend[Lay-1]=Oascend[Lay];
//}
//
//Oascend[Lay-1]=Ohelp*(1.00/(q+1));
//
//UniHelp[0]=Uni[Lay-1];
//IsoHelp[0]=Isometry[Lay-1];
//
//
//int Ihelp=Lay-2;
//Ascen(Oascend, Uni, Rho,Isometry,Ihelp);
//Checker[0](2)=abs(trace(Oascend[Lay-1]*Rho[Lay-1])+trace(Oascend[Lay-1]*Rho[Lay-1]))/2.00;
//
//OptIso(Oascend,Isometry, Uni, Rho,Env, Lay-1);
//OptUni(Oascend,Isometry, Uni, Rho,EnvU, Lay-1);
//
//RhoLanczos(Oascend, Uni, Rho,Isometry);
//
//
// Ihelp=Lay-2;
//Ascen(Oascend, Uni, Rho,Isometry,Ihelp);
//Checker[0](3)=abs(trace(Oascend[Lay-1]*Rho[Lay-1])+trace(Oascend[Lay-1]*Rho[Lay-1]))/2.00;
//if(q==0)
//cout<<setprecision(15)<<"EnergyT="<<Checker[0](2)<<endl;
//if(Checker[0](3)<Checker[0](2)){
//Uni[Lay-1]=UniHelp[0];
//Isometry[Lay-1]=IsoHelp[0];
//RhoLanczos(Oascend, Uni, Rho,Isometry);
//}
//
//}
//
//Ihelp=Lay-2;
//Ascen(Oascend, Uni, Rho,Isometry,Ihelp);
//RhoLanczos(Oascend, Uni, Rho,Isometry);
//
//cout<<"EnergyT="<<abs(trace(Oascend[Lay-1]*Rho[Lay-1])+trace(Oascend[Lay-1]*Rho[Lay-1]))/2.00<<endl;
//
//for(int i=Lay;i>0;i--)
//Dscen(Oascend, Uni, Rho,Isometry,i);
//
//delete [] UniHelp;
//delete [] Checker;
//delete [] IsoHelp;
//}
//
//void OptTopFinite(cx_mat*Oascend,cx_mat*Isometry,cx_mat*Uni,cx_mat*Rho,double &GapEnergy){
//
//cx_mat h;
//cx_mat hB;
//
//cx_mat help;
//cx_mat eigvec;
//vec eigval;
//
//
//hB.zeros(Xi[Lay]*Xi[Lay],Xi[Lay]*Xi[Lay]);
//for(int n=0;n<Xi[Lay];n++)
//for(int m=0;m<Xi[Lay];m++)
//for(int j=0;j<Xi[Lay];j++)
//for(int i=0;i<Xi[Lay];i++)
//hB(i*Xi[Lay]+j,m*Xi[Lay]+n)=Oascend[Lay](j*Xi[Lay]+i,n*Xi[Lay]+m);
//
//h=Oascend[Lay]+hB;
//help=trans(Isometry[Lay])*h*Isometry[Lay];
//
//cout<<"TopEnergy="<<abs(help(0,0)/2.00)<<endl;
//
//eig_sym(eigval, eigvec, h);
//cout<<"TopEnergy="<<abs((1/2.0)*eigval(0))<<endl;
//
//cout<<"Gap="<<pow(3,Lay)*abs(abs(eigval(1))-abs(eigval(0)))<<endl;
//GapEnergy=pow(3,Lay)*abs(abs(eigval(1))-abs(eigval(0)));
//for(int i=0;i<Xi[Lay];i++)
//        for(int j=0;j<Xi[Lay];j++)
//        Isometry[Lay](i*Xi[Lay]+j,0)=eigvec(Xi[Lay]*i+j,0);
//
//help=trans(Isometry[Lay])*h*Isometry[Lay];
////cout<<"TopEnergy="<<help(0,0)/2.00<<endl;
//
//                   for(int n=0;n<Xi[Lay];n++)
//                        for(int m=0;m<Xi[Lay];m++)
//                            for(int j=0;j<Xi[Lay];j++)
//                                for(int i=0;i<Xi[Lay];i++)
//Rho[Lay](i*Xi[Lay]+j,m*Xi[Lay]+n)=Isometry[Lay](i*Xi[Lay]+j,0)*conj(Isometry[Lay](m*Xi[Lay]+n,0));
//
//
//help.zeros(Xi[Lay]*Xi[Lay],Xi[Lay]*Xi[Lay]);
//                    for(int n=0;n<Xi[Lay];n++)
//                        for(int m=0;m<Xi[Lay];m++)
//                            for(int j=0;j<Xi[Lay];j++)
//                                for(int i=0;i<Xi[Lay];i++)
//help(j*Xi[Lay]+i,n*Xi[Lay]+m)=Isometry[Lay](i*Xi[Lay]+j,0)*conj(Isometry[Lay](m*Xi[Lay]+n,0));
//
//Rho[Lay]=(help+Rho[Lay])*(0.5);
//
////cout<<"EnergyTest="<<abs(trace(Oascend[Lay]*Rho[Lay])+trace(Oascend[Lay]*Rho[Lay]))/2.00<<endl;
//
//for(int i=Lay;i>0;i--)
//Dscen(Oascend, Uni, Rho,Isometry,i);
//
//
//}
//
//
