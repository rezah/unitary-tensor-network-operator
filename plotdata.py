from numpy import loadtxt
import matplotlib.pyplot as plt



def plot(Count_list1, variance_list1,Count_list2, variance_list2,Count_list3, variance_list3,Count_list4, variance_list4,q,L,hz_list, E_time1, time1,E_time2, time2,E_time3, time3,E_time4, time4):
 file = open("Data/varianceR_"  + str(q) +".txt", "w")
 for index in range(len(variance_list1)):
     file.write(str(Count_list1[index]) + " " + str(variance_list1[index]) + "\n")
 file.close()
 R=loadtxt("Data/varianceR_"  + str(q) +".txt")

 X=R[:,0]
 Y=R[:,1]



 file = open("Data/varianceM_"  + str(q) +".txt", "w")
 for index in range(len(variance_list2)):
     file.write(str(Count_list2[index]) + "    " + str(variance_list2[index]) + "\n")
 file.close()
 R1=loadtxt("Data/varianceM_"  + str(q) +".txt")

 file = open("Data/hz_list_"  + str(q) +".txt", "w")
 for index in range(len(hz_list)):
     file.write(str(hz_list[index]) + "\n")
 file.close()

 X1=R1[:,0]
 Y1=R1[:,1]

 file = open("Data/varianceM1_"  + str(q) +".txt", "w")
 for index in range(len(variance_list3)):
     file.write(str(Count_list3[index]) + "    " + str(variance_list3[index]) + "\n")
 file.close()
 R2=loadtxt("Data/varianceM1_"  + str(q) +".txt")

 
 X2=R2[:,0]
 Y2=R2[:,1]

 file = open("Data/varianceM2_"  + str(q) +".txt", "w")
 for index in range(len(variance_list4)):
     file.write(str(Count_list4[index]) + "    " + str(variance_list4[index]) + "\n")
 file.close()
 R3=loadtxt("Data/varianceM2_"  + str(q) +".txt")

 
 X3=R3[:,0]
 Y3=R3[:,1]





 plt.plot( X, Y,'bs',label='Regular-architecture-'+str(L) )
 plt.plot( X1, Y1,'g.',label='MERA-'+str(L) )
 plt.plot( X2, Y2,'r>',label='B_MERA_Origin-'+str(L) )
 plt.plot( X3, Y3,'m*',label='T_MERA_Origin-'+str(L) )

 plt.xlabel('Iterations', fontsize=20)
 plt.ylabel(r'$\frac{\bar{\sigma}}{2^{N}}$', fontsize=20)
 plt.legend(loc='upper right')
 plt.savefig('Figs/Convergance'+ str(q) + '.pdf')
 plt.clf()
 
 plt.plot( time2, E_time2,'bs',label='Regular-architecture-'+str(L) )
 plt.plot( time3, E_time3,'g.',label='MERA-'+str(L) )
 plt.plot( time1, E_time1,'r>',label='B_MERA_Origin-'+str(L) )
 plt.plot( time4, E_time4,'m*',label='T_MERA_Origin-'+str(L) )
 plt.xlabel('time', fontsize=20)
 plt.ylabel(r'$\frac{\bar{\sigma}}{2^{N}}$', fontsize=20)
 plt.legend(loc='upper right')
 plt.savefig('Figs/time'+ str(q) + '.pdf')
 plt.clf()


def Store(variance_final1, variance_final2,variance_final3,variance_final4,q,L, file, Count_final1, Count_final2, Count_final3, Count_final4,time1, time2, time3, time4):
 Length=len(variance_final1)-1
 Length2=len(variance_final2)-1
 Length3=len(variance_final3)-1
 Length4=len(variance_final4)-1

 Length11=len(Count_final1)-1
 Length22=len(Count_final2)-1
 Length33=len(Count_final3)-1
 Length44=len(Count_final4)-1
 
 var1=sum(variance_final1)/len(variance_final1)
 var2=sum(variance_final2)/len(variance_final2)
 var3=sum(variance_final3)/len(variance_final3)
 var4=sum(variance_final4)/len(variance_final4)


 
 file.write(str(q)  + " " + str(variance_final1[Length]) + " " + str(variance_final2[Length2]) +  " "+str(variance_final3[Length3])+" "+str(variance_final4[Length4])+" "+ str(var1)+" "+ str(var2) +" "+ str(var3)+" "+str(var4)+" "+ str(Count_final1[Length11])+ " "+ str(Count_final2[Length22])+" "+str(Count_final3[Length33])+"  "+str(Count_final4[Length44])+" "+str(time1[len(time1)-1]) + " "+ str(time2[len(time2)-1])+ " "+ str(time3[len(time3)-1])  +" "+str(time4[len(time4)-1])+  "\n")
 file.flush()
