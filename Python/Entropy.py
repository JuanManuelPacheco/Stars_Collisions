import numpy as np
import sys

folder = sys.argv[1]
fin = sys.argv[2]

# Path where the ascii files of the snapshots are located.
data_path = '/home/juanmanuel/Documents/Master/Thesis/Simulations'

l0=69599000000.000000
m0=1.9891000000000000E+033
d0=m0/(l0*l0*l0)
p0=11253326797072562.
a=7.5646e-15

if folder == 'Renzo':
    lrtab1,mtab1,htab1,hetab1,lTtab1,lPtab1 = np.genfromtxt('{}/PARSEC/star1_RENZO_TAMS_4.15MYR'.format(data_path),
                                                        dtype='float',usecols=(2,1,6,7,3,5) ,skip_header=6,
                                                        comments="#", unpack=True)
else:
    lrtab1,mtab1,htab1,hetab1,lTtab1,lPtab1 = np.genfromtxt('{}/PARSEC/star1_PARSEC'.format(data_path),
                                                        dtype='float',usecols=(2,1,6,7,3,5) ,skip_header=6,
                                                        comments="#", unpack=True)

x,y,z,m,d,e,mu,t = np.genfromtxt('{}/{}/out0000.sph.ascii'.format(data_path,folder),dtype='float',
                                 usecols=(0,1,2,6,7,8,9,12) ,skip_header=12,comments="#", unpack=True)

xf,yf,zf,mf,df,ef,muf,tf = np.genfromtxt('{}/{}/out{}.sph.ascii'.format(data_path,folder,fin),
                                         dtype='float',usecols=(0,1,2,6,7,8,9,12),
                                         skip_header=12,comments="#", unpack=True)

x=x-x[np.argmax(d)]
y=y-y[np.argmax(d)]
z=z-z[np.argmax(d)]

r=np.sqrt(x*x+y*y+z*z)

xf=xf-xf[np.argmax(df)]
yf=yf-yf[np.argmax(df)]
zf=zf-zf[np.argmax(df)]

rf=np.sqrt(xf*xf+yf*yf+zf*zf)

h1,he1 = np.genfromtxt('{}/{}/init_ab_frac1.dat'.format(data_path,folder),dtype='float',
                       usecols=(0,1) ,skip_header=1,comments="#", unpack=True)

h2,he2 = np.genfromtxt('{}/{}/init_ab_frac2.dat'.format(data_path,folder),dtype='float',
                       usecols=(0,1) ,skip_header=1,comments="#", unpack=True)

htot=np.concatenate((h1,h2),axis=None)
hetot=np.concatenate((he1,he2),axis=None)

tens1=lrtab1*0.+10.

rtab1=np.power(tens1,lrtab1)
Ptab1=np.power(tens1,lPtab1)
Ttab1=np.power(tens1,lTtab1)

Ptab1=np.flip(Ptab1)
Ttab1=np.flip(Ttab1)
rtab1=np.flip(rtab1)
mtab1=np.flip(mtab1)

tr=d*0.

if folder == 'Halfres':
    tr[0:449971]=1.

else:
    tr[0:799021]=1.

t=t/4.
tf=tf/4.

d=d/d0
df=df/d0

p=5./3.*d*e*p0
pf=5./3.*df*ef*p0

radp=a*t*t*t*t/3.
radpf=a*tf*tf*tf*tf/3.

ptot=p+radp
pftot=pf+radpf

A=ptot/(np.power(d,5./3.))
Af=pftot/(np.power(df,5./3.))

msort=[x for _,x in sorted(zip(A,m))]
trAsort=[x for _,x in sorted(zip(A,tr))]
hsort=[x for _,x in sorted(zip(A,htot))]
hesort=[x for _,x in sorted(zip(A,hetot))]
Asort=[x for _,x in sorted(zip(A,A))]

trAsort=np.array(trAsort)
hsort=np.array(hsort)
hesort=np.array(hesort)
msort=np.array(msort)
Asort=np.array(Asort)

mcpart=[]
mcum=0.
for i in np.arange(len(msort)):
    mcum=mcum+msort[i]
    mcpart.append(mcum)
mcpart=np.array(mcpart)

mfsort=[x for _,x in sorted(zip(rf,mf))]
Afsort=[x for _,x in sorted(zip(rf,Af))]
trfsort=[x for _,x in sorted(zip(rf,tr))]
hfsort=[x for _,x in sorted(zip(rf,htot))]
hefsort=[x for _,x in sorted(zip(rf,hetot))]
trfsort=np.array(trfsort)
mfsort=np.array(mfsort)
Afsort=np.array(Afsort)
hfsort=np.array(hfsort)
hefsort=np.array(hefsort)

mcfpart=[]
mcum=0.
for i in np.arange(len(mfsort)):
    mcum=mcum+mfsort[i]
    mcfpart.append(mcum)
mcfpart=np.array(mcfpart)

f_bin=open('{}/{}/Entropy.dat'.format(data_path,folder),'w')
f_bin.write('mcfpart Afsort trfsort mcpart Asort trAsort\n')
for i in np.arange(len(mcfpart)):
        f_bin.write(str(mcfpart[i])+' '+str(Afsort[i])+' '+str(trfsort[i])+' '+str(mcpart[i])+' '+str(Asort[i])+' '+str(trAsort[i])+'\n')
f_bin.close()
