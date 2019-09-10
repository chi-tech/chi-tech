import numpy as np
import matplotlib.pyplot as plt
import importlib

execfile('../../ZLFFI00.py')


dz=20
data_avg1 = np.zeros((2,dz))
data_avg2 = np.zeros((2,dz))
err_avg12 = np.zeros((2,dz))
di = 500/dz
for i in range(0,dz):
    data_avg1[0,i] = np.average(data1[(i*di):((i+1)*di),3])
    data_avg1[1,i] = np.average(data1[(i*di):((i+1)*di),4])

    data_avg2[0,i] = np.average(data2[(i*di):((i+1)*di),3])
    data_avg2[1,i] = np.average(data2[(i*di):((i+1)*di),4])

    err_avg12[0,i] = data_avg2[0,i]
    err_avg12[1,i] = -data_avg1[1,i] + data_avg2[1,i]

plt.figure(1)
plt.plot(data0[:,3],data0[:,4],label="Flux_g0")
# plt.plot(data1[:,3],data1[:,4],label="Flux_g0_m0")
# plt.plot(data2[:,3],data2[:,4],label="Flux_g0")
# plt.plot(data_avg1[0,:],data_avg1[1,:])
# plt.plot(data_avg2[0,:],data_avg2[1,:])
plt.plot(err_avg12[0,:],err_avg12[1,:],"kd")
plt.legend()
plt.grid(which='major')
plt.show()
