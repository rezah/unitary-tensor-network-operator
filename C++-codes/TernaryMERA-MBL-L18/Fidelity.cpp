//#include"Header.h"
//void GSF(cx_mat * RhoFidel,cx_mat * RhoFidelp,ofstream & Fidelity,ofstream & Distance,ofstream&Fidelity2,double&J)
//{
//cx_mat RhoFidelpHelp;
//cx_mat DisHelp;
//cx_mat SqrtRho;
//cx_mat FdelityOperator;
//cx_mat eigvec;
//vec eigval;
//for (int i=0;i<Lay ;i++ ){
//eig_sym(eigval, eigvec, RhoFidel[i]);
//eigval=abs(eigval);
//for(int j=0;j<Xi[i]*Xi[i];j++){
//if(eigval(j)<1.0e-15)
//eigval(j)=0;
//}
//SqrtRho=eigvec*sqrt(diagmat(eigval))*trans(eigvec);
//SqrtRho=SqrtRho*(1.00/(sqrt(trace(SqrtRho*SqrtRho))));
//FdelityOperator=SqrtRho*RhoFidelp[i]*SqrtRho;
//eig_sym(eigval, eigvec, FdelityOperator);
//eigval=abs(eigval);
//for(int j=0;j<Xi[i]*Xi[i];j++){
//if(eigval(j)<1.0e-15)
//eigval(j)=0;
//}
//Fidelity<<"  "<<setprecision(15)<<abs(trace(sqrt(diagmat(eigval))));
//DisHelp=RhoFidelp[i]-RhoFidel[i];
//eig_sym(eigval, eigvec, DisHelp);
//eigval=abs(eigval);
//Distance<<"  "<<setprecision(15)<<(1.00-(0.5*trace(diagmat(eigval))));
//Fidelity2<<"  "<<setprecision(15)<<
//abs(abs(trace(RhoFidel[i]*RhoFidelp[i]))*(1.00/(sqrt(  trace(RhoFidel[i]*RhoFidel[i])*trace(RhoFidelp[i]*RhoFidelp[i])    ))));
//}
//Fidelity<<endl;
//Fidelity2<<endl;
//Distance<<endl;
//}
//void GSFMORE(cx_mat * RhoFidel,cx_mat * RhoFidelp,ofstream & Fidelity,ofstream & Distance,ofstream&Fidelity2,double&J)
//{
//cx_mat RhoFidelpHelp;
//cx_mat DisHelp;
//cx_mat SqrtRho;
//cx_mat FdelityOperator;
//cx_mat eigvec;
//vec eigval;
//for (int i=0;i<Lay ;i++ ){
//eig_sym(eigval, eigvec, RhoFidel[i]);
//eigval=abs(eigval);
//for(int j=0;j<Xi[i]*Xi[i];j++){
//if(eigval(j)<1.0e-15)
//eigval(j)=0;
//}
//SqrtRho=eigvec*sqrt(diagmat(eigval))*trans(eigvec);
//SqrtRho=SqrtRho*(1.00/(sqrt(trace(SqrtRho*SqrtRho))));
//FdelityOperator=SqrtRho*RhoFidelp[i]*SqrtRho;
//eig_sym(eigval, eigvec, FdelityOperator);
//eigval=abs(eigval);
//for(int j=0;j<Xi[i]*Xi[i];j++){
//if(eigval(j)<1.0e-15)
//eigval(j)=0;
//}
//Fidelity<<"  "<<setprecision(15)<<abs(trace(sqrt(diagmat(eigval))));
//DisHelp=RhoFidelp[i]-RhoFidel[i];
//eig_sym(eigval, eigvec, DisHelp);
//eigval=abs(eigval);
//Distance<<"  "<<setprecision(15)<<(1.00-(0.5*trace(diagmat(eigval))));
//Fidelity2<<"  "<<setprecision(15)<<
//abs(abs(trace(RhoFidel[i]*RhoFidelp[i]))*(1.00/(sqrt(  trace(RhoFidel[i]*RhoFidel[i])*trace(RhoFidelp[i]*RhoFidelp[i])    ))));
//}
//Fidelity<<endl;
//Fidelity2<<endl;
//Distance<<endl;
//}
//
