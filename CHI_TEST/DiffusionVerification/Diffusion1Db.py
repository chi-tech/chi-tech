import numpy as np
import matplotlib.pyplot as plt

data=np.zeros([25,5])
data[0,0] = 0
data[0,1] = 0
data[0,2] = 0
data[0,3] = 0
data[0,4] = 1.5
data[1,0] = 0.0416667
data[1,1] = 0
data[1,2] = 0
data[1,3] = 0.0416667
data[1,4] = 1.4991
data[2,0] = 0.0833333
data[2,1] = 0
data[2,2] = 0
data[2,3] = 0.0833333
data[2,4] = 1.49647
data[3,0] = 0.125
data[3,1] = 0
data[3,2] = 0
data[3,3] = 0.125
data[3,4] = 1.4921
data[4,0] = 0.166667
data[4,1] = 0
data[4,2] = 0
data[4,3] = 0.166667
data[4,4] = 1.486
data[5,0] = 0.208333
data[5,1] = 0
data[5,2] = 0
data[5,3] = 0.208333
data[5,4] = 1.47817
data[6,0] = 0.25
data[6,1] = 0
data[6,2] = 0
data[6,3] = 0.25
data[6,4] = 1.4686
data[7,0] = 0.291667
data[7,1] = 0
data[7,2] = 0
data[7,3] = 0.291667
data[7,4] = 1.4573
data[8,0] = 0.333333
data[8,1] = 0
data[8,2] = 0
data[8,3] = 0.333333
data[8,4] = 1.44427
data[9,0] = 0.375
data[9,1] = 0
data[9,2] = 0
data[9,3] = 0.375
data[9,4] = 1.4295
data[10,0] = 0.416667
data[10,1] = 0
data[10,2] = 0
data[10,3] = 0.416667
data[10,4] = 1.413
data[11,0] = 0.458333
data[11,1] = 0
data[11,2] = 0
data[11,3] = 0.458333
data[11,4] = 1.39477
data[12,0] = 0.5
data[12,1] = 0
data[12,2] = 0
data[12,3] = 0.5
data[12,4] = 1.3748
data[13,0] = 0.541667
data[13,1] = 0
data[13,2] = 0
data[13,3] = 0.541667
data[13,4] = 1.3531
data[14,0] = 0.583333
data[14,1] = 0
data[14,2] = 0
data[14,3] = 0.583333
data[14,4] = 1.32967
data[15,0] = 0.625
data[15,1] = 0
data[15,2] = 0
data[15,3] = 0.625
data[15,4] = 1.3045
data[16,0] = 0.666667
data[16,1] = 0
data[16,2] = 0
data[16,3] = 0.666667
data[16,4] = 1.2776
data[17,0] = 0.708333
data[17,1] = 0
data[17,2] = 0
data[17,3] = 0.708333
data[17,4] = 1.24897
data[18,0] = 0.75
data[18,1] = 0
data[18,2] = 0
data[18,3] = 0.75
data[18,4] = 1.2186
data[19,0] = 0.791667
data[19,1] = 0
data[19,2] = 0
data[19,3] = 0.791667
data[19,4] = 1.1865
data[20,0] = 0.833333
data[20,1] = 0
data[20,2] = 0
data[20,3] = 0.833333
data[20,4] = 1.15267
data[21,0] = 0.875
data[21,1] = 0
data[21,2] = 0
data[21,3] = 0.875
data[21,4] = 1.1171
data[22,0] = 0.916667
data[22,1] = 0
data[22,2] = 0
data[22,3] = 0.916667
data[22,4] = 1.0798
data[23,0] = 0.958333
data[23,1] = 0
data[23,2] = 0
data[23,3] = 0.958333
data[23,4] = 1.04077
data[24,0] = 1
data[24,1] = 0
data[24,2] = 0
data[24,3] = 1
data[24,4] = 1
done=True


x = np.linspace(0,1,25)
y = np.zeros(25)
for i in range(0,25):
    y[i] = -0.5*x[i]*x[i] + 3/2


plt.figure(1)
plt.plot(data[:,3],data[:,4],'kD',label='CFEM 1D',fillstyle='none')
plt.plot(x,y,'k-',label='Analytical solution')

plt.legend()
plt.xlabel('x')
plt.ylabel(r'$\phi$')

plt.show()

