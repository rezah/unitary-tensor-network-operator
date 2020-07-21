#include"Header.h"
void UniInitial(cx_mat*Uni0,cx_mat*Uni1,cx_mat *Isometry0,cx_mat *Isometry1)
{
cx_mat d1;
cx_mat d2;
cx_mat C1;
cx_mat C2;
cx_mat V1;
vec s1;


//arma_rng::set_seed_random();

for(int i=0;i<int(N0/3);i++){
C1=randu<cx_mat>(Xi0[1]*Xi0[2]*Xi0[0],Xi0[1]*Xi0[2]*Xi0[0]);
svd_econ(d1, s1, V1, C1,'l',"std");
Isometry0[i]=d1;
Isometry0[i].eye(Xi0[2]*Xi0[2]*Xi0[0],Xi0[2]*Xi0[2]*Xi0[0]);
}

for(int i=0;i<int(N0/3)-1;i++){
C1=randu<cx_mat>(Xi0[2]*Xi0[2],Xi0[2]*Xi0[2]);
svd_econ(V1, s1, d1, C1,'b',"std");
Uni0[i]=V1;
Uni0[i].eye(Xi0[2]*Xi0[2],Xi0[2]*Xi0[2]);
}



C1=randn<cx_mat>(Xi1[2]*Xi1[3],Xi1[2]*Xi1[3]);
svd_econ(d1, s1, V1, C1,'l',"std");
Uni1[0]=d1;
Uni1[0].eye(Xi1[2]*Xi1[3],Xi1[2]*Xi1[3]);

C1=randn<cx_mat>(Xi1[0]*Xi1[1]*Xi1[2],Xi1[0]*Xi1[1]*Xi1[2]);
svd_econ(d1, s1, V1, C1,'l',"std");
Isometry1[0]=d1;
Isometry1[0].eye(Xi1[0]*Xi1[1]*Xi1[2],Xi1[0]*Xi1[1]*Xi1[2]);


C1=randn<cx_mat>(Xi1[3]*Xi1[4]*Xi1[5],Xi1[3]*Xi1[4]*Xi1[5]);
svd_econ(d1, s1, V1, C1,'l',"std");
Isometry1[1]=d1;
Isometry1[1].eye(Xi1[3]*Xi1[4]*Xi1[5],Xi1[3]*Xi1[4]*Xi1[5]);






}


