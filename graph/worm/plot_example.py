import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-whitegrid')
Avr1=[12,11,7,7,6,5]
Var1=[0.5,0.4,0.3,1,0.3,0.5]
Avr2=[10,8,5,4,3,3]
Var2=[0.4,0.3,0.4,0.6,0.3,0.5]
Time=np.arange(1,7,1)
plt.errorbar(Time,Avr1,yerr=Var1,fmt='.k')
plt.errorbar(Time,Avr2,yerr=Var2,capsize=2)
plt.xlabel('moth')
plt.ylabel('ame/cm') 
plt.show()