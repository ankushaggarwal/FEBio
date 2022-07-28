import pickle
import struct

ff = open('splines.p','rb')
I1 = pickle.load(ff)
I4 = pickle.load(ff)
mean_sp = pickle.load(ff)
ff.close()

ff = open('example.dat','wb')
degrees = mean_sp.degrees
k = mean_sp.get_knots()
c = mean_sp.get_coeffs()

for j in range(2):
    ff.write(struct.pack('=ii',degrees[j],len(k[j])))
    for i in range(len(k[j])):
        ff.write(struct.pack('=d',k[j][i]))

ff.write(struct.pack('=i',len(c)))
for i in range(len(c)):
    ff.write(struct.pack('=d',c[i]))
ff.close()

