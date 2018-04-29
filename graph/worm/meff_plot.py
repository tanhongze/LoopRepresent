from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
from data import meff_b as mb
if __name__ == '__main__':
    fig = plt.figure()
    plt.plot(np.arange(0.3,-0.59,-0.02),mb)
    plt.show()