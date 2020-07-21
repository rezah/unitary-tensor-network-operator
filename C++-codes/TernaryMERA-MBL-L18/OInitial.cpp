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
Oascend1[0]=kron(kron(kron(PX,PX),PI),kron(kron(PI,PI),PI))+kron(kron(kron(PY,PY),PI),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PZ,PZ),PI),kron(kron(PI,PI),PI))+Randomvec(0)*kron(kron(kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[0]=Oascend1[0]+kron(kron(kron(PI,PX),PX),kron(kron(PI,PI),PI))+kron(kron(kron(PI,PY),PY),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PI,PZ),PZ),kron(kron(PI,PI),PI))+Randomvec(1)*kron(kron(kron(PI,PZ),PI),kron(kron(PI,PI),PI));

Oascend1[0]=Oascend1[0]+kron(kron(kron(PI,PI),PX),kron(kron(PX,PI),PI))+kron(kron(kron(PI,PI),PY),kron(kron(PY,PI),PI))+Jz*kron(kron(kron(PI,PI),PZ),kron(kron(PZ,PI),PI))+Randomvec(2)*kron(kron(kron(PI,PI),PZ),kron(kron(PI,PI),PI));

///////////
Oascend1[1]=kron(kron(kron(PX,PX),PI),kron(kron(PI,PI),PI))+kron(kron(kron(PY,PY),PI),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PZ,PZ),PI),kron(kron(PI,PI),PI))+Randomvec(3)*kron(kron(kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[1]=Oascend1[1]+kron(kron(kron(PI,PX),PX),kron(kron(PI,PI),PI))+kron(kron(kron(PI,PY),PY),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PI,PZ),PZ),kron(kron(PI,PI),PI))+Randomvec(4)*kron(kron(kron(PI,PZ),PI),kron(kron(PI,PI),PI));

Oascend1[1]=Oascend1[1]+kron(kron(kron(PI,PI),PX),kron(kron(PX,PI),PI))+kron(kron(kron(PI,PI),PY),kron(kron(PY,PI),PI))+Jz*kron(kron(kron(PI,PI),PZ),kron(kron(PZ,PI),PI))+Randomvec(5)*kron(kron(kron(PI,PI),PZ),kron(kron(PI,PI),PI));
//////////

Oascend1[2]=kron(kron(kron(PX,PX),PI),kron(kron(PI,PI),PI))+kron(kron(kron(PY,PY),PI),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PZ,PZ),PI),kron(kron(PI,PI),PI))+Randomvec(6)*kron(kron(kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[2]=Oascend1[2]+kron(kron(kron(PI,PX),PX),kron(kron(PI,PI),PI))+kron(kron(kron(PI,PY),PY),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PI,PZ),PZ),kron(kron(PI,PI),PI))+Randomvec(7)*kron(kron(kron(PI,PZ),PI),kron(kron(PI,PI),PI));

Oascend1[2]=Oascend1[2]+kron(kron(kron(PI,PI),PX),kron(kron(PX,PI),PI))+kron(kron(kron(PI,PI),PY),kron(kron(PY,PI),PI))+Jz*kron(kron(kron(PI,PI),PZ),kron(kron(PZ,PI),PI))+Randomvec(8)*kron(kron(kron(PI,PI),PZ),kron(kron(PI,PI),PI));
///////////

Oascend1[3]=kron(kron(kron(PX,PX),PI),kron(kron(PI,PI),PI))+kron(kron(kron(PY,PY),PI),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PZ,PZ),PI),kron(kron(PI,PI),PI))+Randomvec(9)*kron(kron(kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[3]=Oascend1[3]+kron(kron(kron(PI,PX),PX),kron(kron(PI,PI),PI))+kron(kron(kron(PI,PY),PY),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PI,PZ),PZ),kron(kron(PI,PI),PI))+Randomvec(10)*kron(kron(kron(PI,PZ),PI),kron(kron(PI,PI),PI));

Oascend1[3]=Oascend1[3]+kron(kron(kron(PI,PI),PX),kron(kron(PX,PI),PI))+kron(kron(kron(PI,PI),PY),kron(kron(PY,PI),PI))+Jz*kron(kron(kron(PI,PI),PZ),kron(kron(PZ,PI),PI))+Randomvec(11)*kron(kron(kron(PI,PI),PZ),kron(kron(PI,PI),PI));
///////

Oascend1[4]=kron(kron(kron(PX,PX),PI),kron(kron(PI,PI),PI))+kron(kron(kron(PY,PY),PI),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PZ,PZ),PI),kron(kron(PI,PI),PI))+Randomvec(12)*kron(kron(kron(PZ,PI),PI),kron(kron(PI,PI),PI));

Oascend1[4]=Oascend1[4]+kron(kron(kron(PI,PX),PX),kron(kron(PI,PI),PI))+kron(kron(kron(PI,PY),PY),kron(kron(PI,PI),PI))+Jz*kron(kron(kron(PI,PZ),PZ),kron(kron(PI,PI),PI))+Randomvec(13)*kron(kron(kron(PI,PZ),PI),kron(kron(PI,PI),PI));


Oascend1[4]=Oascend1[4]+kron(kron(kron(PI,PI),PX),kron(kron(PX,PI),PI))+kron(kron(kron(PI,PI),PY),kron(kron(PY,PI),PI))+Jz*kron(kron(kron(PI,PI),PZ),kron(kron(PZ,PI),PI))+Randomvec(14)*kron(kron(kron(PI,PI),PZ),kron(kron(PI,PI),PI));


Oascend1[4]=Oascend1[4]+kron(kron(kron(PI,PI),PI),kron(kron(PX,PX),PI))+kron(kron(kron(PI,PI),PI),kron(kron(PY,PY),PI))+Jz*kron(kron(kron(PI,PI),PI),kron(kron(PZ,PZ),PI))+Randomvec(15)*kron(kron(kron(PI,PI),PI),kron(kron(PZ,PI),PI))+Randomvec(16)*kron(kron(kron(PI,PI),PI),kron(kron(PI,PZ),PI));


Oascend1[4]=Oascend1[4]+kron(kron(kron(PI,PI),PI),kron(kron(PI,PX),PX))+kron(kron(kron(PI,PI),PI),kron(kron(PI,PY),PY))+Jz*kron(kron(kron(PI,PI),PI),kron(kron(PI,PZ),PZ))+Randomvec(17)*kron(kron(kron(PI,PI),PI),kron(kron(PI,PI),PZ));


cout<<((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
cout<<(2*(N0-1)/16.00)+((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00)<<endl;
traceH2=((2*(N0-1)/16.00)+((Jz*Jz*(N0-1))/16.00)+(trace(trans(Randomvec)*Randomvec)/4.00));

}



