import os, math
import numpy as np
import matplotlib.pyplot as plt

rad2deg = lambda x : [i*180.0/math.pi for i in x]

for filename in os.listdir('.'):
    if filename == 'Exam4_matlab_output.txt':
        time, q1, q2, u1, u2 = np.loadtxt(open(filename),
                                unpack=True,
                                skiprows=5)
        q1 = rad2deg(q1)
        q2 = rad2deg(q2)
        finalQ1 = q1[-1]
        finalQ2 = q2[-1]
        plt.plot(time,q1, label='Q1 Final Angle = {:.4f}$^\circ$'.format(finalQ1))
        plt.plot(time,q2, label='Q2 Final Angle = {:.4f}$^\circ$'.format(finalQ2))
        plt.xlabel('Time (s)')
        plt.ylabel('Angle (degrees)')
        plt.legend()
        plt.title('Restrained Double Pendulum')
        plt.savefig(filename[:-4]+'_figure.png', bbox_inches='tight')
        plt.show()
