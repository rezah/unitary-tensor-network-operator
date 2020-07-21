#include"Header.h"
void OInitial(cx_mat*Oascend0,cx_mat*Oascend1,double&J,double & traceH2)
{
complex<double> img(0,1);
vec eigval;
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
PX=PX*(1.00/2.00);
PY=PY*(1.00/2.00);
PZ=PZ*(1.00/2.00);

/*************************/
vec Randomvec=randu<vec>(int(N0));

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);

  std::uniform_real_distribution<double> distribution (int(-J),int(J));

  std::cout << "some random numbers between -J and J: "<<J<< std::endl;
  for (int i=0; i<int(N0); ++i){
    Randomvec(i)=distribution(generator);
    }
cout<<Randomvec<<endl;
//arma_rng::set_seed_random();
//for(int i=0;i<int(N0);i++)
//cout<<J*(2*randu(1)-1)<<endl;

if(NewGeneration==1)
//Randomvec.save("Storage/Randomvec.vec",arma_binary);

Randomvec.load("Storage/Randomvec.vec",arma_binary);

//   for (int i=0; i<int(N0); ++i){
//     Randomvec(i)=1.0;
//     }

cout<<endl<<Randomvec<<endl;

//for(int i=0;i<int(N0);i++)
//Oascend0[i]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)+Randomvec(i)*kron(PZ,PI);                     //Ising

//cout<<"hi"<<endl;
//for(int i=0;i<int(N0)-1;i++){
//Oascend1[i]=kron(PX,PX)+kron(PY,PY)+kron(PZ,PZ)+Randomvec(i)*kron(PZ,PI);                     //Ising
//if(i==int(N0)-2)
//Oascend1[i]=Randomvec(i+1)*kron(PI,PZ);
//}
double Jz;
Jz=0.0;
Oascend1[0]=kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(0)*kron(PZ,PI);
Oascend1[1]=kron(kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(1)*kron(PZ,PI),kron(PI,PI));

Oascend1[2]=kron(kron(kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(2)*kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[2]=kron(kron(PI,kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(3)*kron(PZ,PI)),kron(kron(PI,PI),PI))+Oascend1[2];

Oascend1[2]=kron(kron(kron(PI,PI),kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(4)*kron(PZ,PI)),kron(PI,PI))+
Oascend1[2];

Oascend1[2]=kron(kron(kron(PI,PI),PI),kron(kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(5)*kron(PZ,PI),PI))+Oascend1[2];

Oascend1[2]=kron(kron(kron(PI,PI),PI),kron(PI,kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(6)*kron(PZ,PI)))+Oascend1[2];


Oascend1[3]=kron(kron(PI,PI),kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(7)*kron(PZ,PI));
Oascend1[4]=kron(PX,PX)+kron(PY,PY)+Jz*kron(PZ,PZ)+Randomvec(8)*kron(PZ,PI)+Randomvec(9)*kron(PI,PZ);

cout<<((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
cout<<(2*(N0-1)/16.00)+((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
traceH2=((2*(N0-1)/16.00)+((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00));

}



