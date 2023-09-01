import numpy as np
import sys

folder = sys.argv[1]

# Path where the ascii files of the snapshots are located.
data_path = '/home/juanmanuel/Documents/Master/Thesis/Simulations'

if folder == 'Renzo':
    lrtab1,mtab1,htab1,hetab1,lTtab1,lPtab1 = np.genfromtxt('{}/PARSEC/star1_RENZO_TAMS_4.15MYR'.format(data_path),
                                                        dtype='float',usecols=(2,1,6,7,3,5) ,skip_header=6,
                                                        comments="#", unpack=True)
else:
    lrtab1,mtab1,htab1,hetab1,lTtab1,lPtab1 = np.genfromtxt('{}/PARSEC/star1_PARSEC'.format(data_path),
                                                        dtype='float',usecols=(2,1,6,7,3,5) ,skip_header=6,
                                                        comments="#", unpack=True)

if folder == '4Myr':
    lrtab2,mtab2,htab2,hetab2 = np.genfromtxt('{}/PARSEC/star2_PARSEC_MS_4.0MYR'.format(data_path),
                                          dtype='float', usecols=(2,1,6,7) ,skip_header=6,
                                          comments="#", unpack=True)

else:
    lrtab2,mtab2,htab2,hetab2 = np.genfromtxt('{}/PARSEC/star2_PARSEC'.format(data_path),
                                          dtype='float', usecols=(2,1,6,7) ,skip_header=6,
                                          comments="#", unpack=True)

tens1=lrtab1*0.+10.
tens2=lrtab2*0.+10.


rtab1=np.power(tens1,lrtab1)
rtab2=np.power(tens2,lrtab2)

rtab1=np.flip(rtab1)
htab1=np.flip(htab1)
hetab1=np.flip(hetab1)
rtab2=np.flip(rtab2)
htab2=np.flip(htab2)
hetab2=np.flip(hetab2)
mtab1=np.flip(mtab1)
mtab2=np.flip(mtab2)


x,y,z,m,d,e,mu,t = np.genfromtxt('{}/{}/out0000.sph.ascii'.format(data_path,folder),dtype='float',
                                 usecols=(0,1,2,6,7,8,9,12) ,skip_header=12,comments="#", unpack=True)

x=x-x[np.argmax(d)]
y=y-y[np.argmax(d)]
z=z-z[np.argmax(d)]

r=np.sqrt(x*x+y*y+z*z)

index = np.argsort(r)

rsort = r[index]
msort = m[index]

h1=r*0.
he1=r*0.

mr=0.
mcum=[]

indm=np.arange(len(r))*0

mcumpart=r*0.

for n in np.arange(0,len(rtab1)):
    
    if(n>0):
        ind1=((r<rtab1[n]) & (r>=rtab1[n-1]))
    else:
        ind1=((r<rtab1[n]))
    
    mr=mr+np.sum(m[ind1])
    mcum.append(mr)
    mcumpart[ind1]=mr

h2=r*0.
he2=r*0.

mr2=0.
mcum2=[]

indm2=np.arange(len(r))*0

mcumpart2=r*0.

for n in np.arange(0,len(rtab2)):
    
    if(n>0):
        ind2=((r<rtab2[n]) & (r>=rtab2[n-1]))
    else:
        ind2=((r<rtab2[n]))
    
    mr2=mr2+np.sum(m[ind2])
    mcum2.append(mr2)
    mcumpart2[ind2]=mr2

fhpart=r*0.
fhepart=r*0.

dm1=mtab1*0.
dm1[0]=mtab1[0]

for i in np.arange(1,len(mtab1)):
    dm1[i]=mtab1[i]-mtab1[i-1]

fhpart2=r*0.
fhepart2=r*0.

dm2=mtab2*0.
dm2[0]=mtab2[0]

for i in np.arange(1,len(mtab2)):
    dm2[i]=mtab2[i]-mtab2[i-1]

for i in np.arange(0,len(mcum)):
    
    if(i>0):
        indcum=((mtab1>mcum[i-1]) & (mtab1<=mcum[i]))
        indpart=((mcumpart>mcum[i-1]) & (mcumpart<=mcum[i]))
    else:
        indcum=((mtab1<=mcum[i]))
        indpart=((mcumpart<=mcum[i]))
    
    mhsh=np.sum(dm1[indcum]*htab1[indcum])
    mhesh=np.sum(dm1[indcum]*hetab1[indcum])
    msh=np.sum(dm1[indcum])    
    
    if(msh!=0.):
        fhsh=mhsh/msh
        fhesh=mhesh/msh
        
    fhpart[indpart]=fhsh
    fhepart[indpart]=fhesh

for i in np.arange(0,len(mcum2)):
    
    if(i>0):
        indcum2=((mtab2>mcum2[i-1]) & (mtab2<=mcum2[i]))
        indpart2=((mcumpart2>mcum2[i-1]) & (mcumpart2<=mcum2[i]))
    else:
        indcum2=((mtab2<=mcum2[i]))
        indpart2=((mcumpart2<=mcum2[i]))
    
    mhsh2=np.sum(dm2[indcum2]*htab2[indcum2])
    mhesh2=np.sum(dm2[indcum2]*hetab2[indcum2])
    msh2=np.sum(dm2[indcum2])    
    
    if(msh2!=0.):
        fhsh2=mhsh2/msh2
        fhesh2=mhesh2/msh2
        
    fhpart2[indpart2]=fhsh2
    fhepart2[indpart2]=fhesh2

f_bin=open('{}/{}/init_ab_frac1.dat'.format(data_path,folder),'w')
f_bin.write('Hfrac1 Hefrac1\n')
for i in np.arange(len(fhpart)):
        f_bin.write(str(fhpart[i])+' '+str(fhepart[i])+' '+str(mcumpart[i])+' '+str(r[i])+'\n')
f_bin.close()

f_bin2=open('{}/{}/init_ab_frac2.dat'.format(data_path,folder),'w')
f_bin2.write('Hfrac2 Hefrac2\n')
for i in np.arange(len(fhpart2)):
        f_bin2.write(str(fhpart2[i])+' '+str(fhepart2[i])+' '+str(mcumpart2[i])+'\n')
f_bin2.close()

