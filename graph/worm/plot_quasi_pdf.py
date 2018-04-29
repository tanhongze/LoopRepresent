from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
from worm_quasi_pdf_611 import pdfs
for keys in pdfs:
    if keys[3]!=0:
        continue
    if keys[4]!=0:
        continue
    if keys[8]!=0:
        continue
    if keys[9]!=1:
        continue
    plt.plot(np.arange(len(pdfs[keys][8])),np.array(pdfs[keys][8])/pdfs[keys][2])
plt.show()