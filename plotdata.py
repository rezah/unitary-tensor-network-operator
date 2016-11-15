from numpy import loadtxt
import matplotlib.pyplot as plt


def plot(Count_list1, variance_list1,Count_list2, variance_list2,q,L,hz_list):
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
 plt.plot( X, Y,'bs',label='Regular-architecture-'+str(L) )
 plt.plot( X1, Y1,'g.',label='MERA-'+str(L) )
 plt.xlabel('Iterations', fontsize=20)
 plt.ylabel(r'$\frac{\bar{\sigma}}{2^{N}}$', fontsize=20)
 plt.legend(loc='upper right')
 plt.savefig('Figs/Convergance'+ str(q) + '.pdf')
 plt.clf()
 
def Store(variance_final1, variance_final2,q,L, file, Count_final1, Count_final2):
 Length=len(variance_final1)-1
 Length1=len(Count_final1)-1
 
 var1=sum(variance_final1)/len(variance_final1)
 var2=sum(variance_final2)/len(variance_final2)
 
 file.write(str(q)  + "  " + str(variance_final1[Length]) + "  " + str(variance_final2[Length]) +  "  "+ str(var1)+"  "+ str(var2) + "  "+ str(Count_final1[Length1])+ "  "+ str(Count_final2[Length1]) +"\n")
 file.flush()
