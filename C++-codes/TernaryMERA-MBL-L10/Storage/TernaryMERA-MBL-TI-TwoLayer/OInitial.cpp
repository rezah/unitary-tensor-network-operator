#include"Header.h"
void OInitial(cx_mat*Oascend0,cx_mat*Oascend1,cx_mat *Uni0,cx_mat*Isometry0 ,double&J,double & traceH2,ofstream &ConfEnergy,
vec * CouplingRandoms)
{
complex<double> img(0,1);
vec eigval;
vec Randomvec=randu<vec>(int(Nreal));

cx_mat PZ(2,2); PZ.zeros();        PZ(0,0)=1; PZ(1,1)=-1;
cx_mat PI(2,2);  PI.eye(2,2);
cx_mat PX(2,2); PX.zeros();    PX(0,1)=1; PX(1,0)=1;
cx_mat PY(2,2); PY.zeros();    PY(0,1)=-img; PY(1,0)=img;
cx_mat I(3,3);  I.eye(3,3);
cx_mat X(3,3);  X.zeros();           X(0,1)=1; X(1,2)=1; X(2,0)=1;
cx_mat Z(3,3);  Z.zeros();           Z(0,0)=2; Z(1,1)=-1; Z(2,2)=-1;
cx_mat X1(3,3);  X1.zeros();         X1(0,1)=1; X1(1,0)=1; X1(1,2)=1,X1(2,1)=1;
cx_mat Z1(3,3);  Z1.zeros();         Z1(0,0)=1; Z1(1,1)=0; Z1(2,2)=-1;
cx_mat Y1(3,3);  Y1.zeros();         Y1(0,1)=-1; Y1(1,0)=1; Y1(1,2)=-1,Y1(2,1)=1;
/*************************/
PX=PX*(1.00/2.00);
PY=PY*(1.00/2.00);
PZ=PZ*(1.00/2.00);



unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);

  std::uniform_real_distribution<double> distribution (int(-J),int(J));

  std::cout << "some random numbers between -J and J: " << std::endl;
  for (int i=0; i<int(Nreal); ++i){
    Randomvec(i)=distribution(generator);
}
cout<<Randomvec<<endl;
//arma_rng::set_seed_random();
//for(int i=0;i<int(N0);i++)
//cout<<J*(2*randu(1)-1)<<endl;
if(NewGeneration==0)
Randomvec.save("Storage/Randomvec.vec",arma_binary);
Randomvec.load("Storage/Randomvec.vec",arma_binary);

//Randomvec.fill(-20);

cout<<endl<<Randomvec<<endl;

CouplingRandoms[0]=Randomvec;

for (int i=0; i<int(Nreal); ++i){
ConfEnergy<<J<<"  "<<Randomvec(i)<<endl;
}



for(int i=1;i<=int(Nreal);i++){
Oascend1[i]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(i-1)*kron(PZ,PI);
}

Oascend1[6]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(0)*kron(PZ,PI);
Oascend1[1]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(1)*kron(PZ,PI);
Oascend1[2]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(2)*kron(PZ,PI);
Oascend1[3]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(3)*kron(PZ,PI);
Oascend1[4]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(4)*kron(PZ,PI);
Oascend1[5]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)-Randomvec(5)*kron(PZ,PI);



cout<<((Nreal)/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
cout<<(3*(Nreal)/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
traceH2=((3*(Nreal)/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00));


cout<<traceH2<<endl;


complex<double>
Htopowerto=trace(Oascend1[1]*Oascend1[1])*(double(Xi1[3]*Xi1[4])*Xi1[5]*Xi1[0]);
Htopowerto=Htopowerto+trace(Oascend1[2]*Oascend1[2])*(double(Xi1[4]*Xi1[5])*Xi1[0]*Xi1[1]);
Htopowerto=Htopowerto+trace(Oascend1[3]*Oascend1[3])*(double(Xi1[5]*Xi1[0])*Xi1[1]*Xi1[2]);
Htopowerto=Htopowerto+trace(Oascend1[4]*Oascend1[4])*(double(Xi1[0]*Xi1[1])*Xi1[2]*Xi1[3]);
Htopowerto=Htopowerto+trace(Oascend1[5]*Oascend1[5])*(double(Xi1[1]*Xi1[2])*Xi1[3]*Xi1[4]);
Htopowerto=Htopowerto+trace(Oascend1[6]*Oascend1[6])*(double(Xi1[2]*Xi1[3])*Xi1[4]*Xi1[5]);

cout<<real(Htopowerto*(1.00/pow(2,Nreal)))<<endl;
traceH2=real(Htopowerto*(1.00/pow(2,Nreal)));

}



