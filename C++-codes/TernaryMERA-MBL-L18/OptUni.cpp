//#include"Header.h"
//void OptUni(cx_mat*Oascend,cx_mat*Isometry,cx_mat*Uni,cx_mat*Rho,cx_mat*EnvU,int NumLay)
//{
//cx_mat U;
//cx_mat V;
//vec s;
//vec eigvalA;
//cx_mat Identity;
//Identity.eye(Xi[NumLay]*Xi[NumLay],Xi[NumLay]*Xi[NumLay]);
//cx_mat * UniHelp=new cx_mat[1];
//
//double CheckEnergy;
//double CheckEnergyF;
//int Xl=Xi[NumLay];
//
//
//UniHelp[0]=Uni[NumLay];
//
//EnvUF(Oascend,Uni,Rho,Isometry,EnvU,NumLay);
//
//CheckEnergy=abs(trace(EnvU[NumLay]*Uni[NumLay])+trace(EnvU[NumLay]*Uni[NumLay]))/2.00;
//
////cout<<setprecision(15)<<endl<<"EnvUni="<<abs(CheckEnergy)<<endl;
//
//eig_sym(eigvalA, Oascend[NumLay]);
//Oascend[NumLay]=Oascend[NumLay]-(eigvalA(Xl*Xl-1))*Identity;
//
//
//for(int z=0;z<10;z++){
//EnvUF(Oascend,Uni,Rho,Isometry,EnvU,NumLay);
//svd_econ(U, s, V, EnvU[NumLay],'b',"std");
//Uni[NumLay]=-V*trans(U);
//}
//Oascend[NumLay]=Oascend[NumLay]+(eigvalA(Xl*Xl-1))*Identity;
//
//EnvUF(Oascend,Uni,Rho,Isometry,EnvU,NumLay);
//
//CheckEnergyF=abs(trace(EnvU[NumLay]*Uni[NumLay])+trace(EnvU[NumLay]*Uni[NumLay]))/2.00;
////cout<<setprecision(16)<<"EnvUni="<<abs(CheckEnergyF)<<endl<<endl;
//
//if(CheckEnergyF<CheckEnergy){
//Uni[NumLay]=UniHelp[0];
//cout<<"Not being optimized, previous result"<<endl;}
///**************************************************************/
//
//
//delete [] UniHelp;
//}
//void OptIso(cx_mat*Oascend,cx_mat*Isometry,cx_mat*Uni,cx_mat*Rho,cx_mat*Env,int NumLay)
//{
//cx_mat U;
//cx_mat V;
//vec s;
//vec eigvalA;
//cx_mat Identity;
//Identity.eye(Xi[NumLay]*Xi[NumLay],Xi[NumLay]*Xi[NumLay]);
//cx_mat * IsometryHelp=new cx_mat[1];
//double CheckEnergy;
//double CheckEnergyF;
//
////cout<<"here="<<NumLay<<endl;
////cout<<Isometry[NumLay].n_cols<<"  "<<Isometry[NumLay].n_rows<<endl;
//int Xl=Xi[NumLay];
//
//IsometryHelp[0]=Isometry[NumLay];
//
//EnvF(Oascend,Uni,Rho,Isometry,Env,NumLay);
//
//CheckEnergy=abs(trace(Isometry[NumLay]*Env[NumLay]));
//
////cout<<setprecision(16)<<endl<<"Env="<<abs(CheckEnergy)<<endl;
//
//eig_sym(eigvalA, Oascend[NumLay]);
//Oascend[NumLay]=Oascend[NumLay]-(eigvalA(Xl*Xl-1))*Identity;
//
//
//for(int z=0;z<10;z++){
//EnvF(Oascend,Uni,Rho,Isometry,Env,NumLay);
//svd_econ(U, s, V, Env[NumLay],'b',"std");
//Isometry[NumLay]=-V*trans(U);
//}
//Oascend[NumLay]=Oascend[NumLay]+(eigvalA(Xl*Xl-1))*Identity;
//
//EnvF(Oascend,Uni,Rho,Isometry,Env,NumLay);
//CheckEnergyF=abs(trace(Isometry[NumLay]*Env[NumLay]));
//
//
////cout<<setprecision(16)<<"Env="<<abs(CheckEnergyF)<<endl<<endl;
////cout<<"here"<<endl;
////cout<<Isometry[NumLay].n_cols<<"  "<<Isometry[NumLay].n_rows<<endl;
//
//if(CheckEnergyF<CheckEnergy){
//Isometry[NumLay]=IsometryHelp[0];
//cout<<"Not being optimized, previous result"<<endl;}
//
//delete [] IsometryHelp;
//}
