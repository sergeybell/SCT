import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = np.loadtxt('WIEN.dat')
#x = data[:,0]
#y = data[:,1]
#plt.plot(x, y)
#plt.show()

df = pd.DataFrame(data)
df.columns = ['EB', 'CONS(MU)', 'CONS(NBAR(1))', 'CONS(NBAR(2))', 'CONS(NBAR(3))']
#df.plot(x='EB', y='CONS(MU)')
#df.plot(x='EB', y='CONS(NBAR(1))')
#df.plot(x='EB', y='CONS(NBAR(2))')
#df.plot(x='EB', y='CONS(NBAR(3))')

fig, ax = plt.subplots(4,1,sharex=True)
t = df['EB']
ax[0].plot(t, df['CONS(MU)']); ax[0].set_ylabel('CONS(MU)')
ax[1].plot(t, df['CONS(NBAR(1))']); ax[1].set_ylabel('CONS(NBAR(1))')
ax[2].plot(t, df['CONS(NBAR(2))']); ax[2].set_ylabel('CONS(NBAR(2))')
ax[3].plot(t, df['CONS(NBAR(3))']); ax[3].set_ylabel('CONS(NBAR(3))')
ax[3].set_xlabel('EB')
plt.show()

#print(df.tail())
#print(df.info())
#df.to_csv('myarray.csv', index=False, header=False)

#df.plot(x="Name", y="Age", kind="bar")


#xpoints = np.array([1, 2, 6, 8])
#ypoints = np.array([3, 8, 1, 10])

#plt.plot(xpoints, ypoints)
#plt.show()
