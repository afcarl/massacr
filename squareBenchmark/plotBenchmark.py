import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
plt.rcParams['contour.negative_linestyle'] = 'solid'


t = np.loadtxt('t.txt',delimiter='\n')
x0 = np.loadtxt('x.txt',delimiter='\n')
x0=x0/np.max(x0)
y0 = np.loadtxt('y.txt',delimiter='\n')
y0=y0/np.max(y0)
x=x0
y=y0
print y.shape
bits = 81
x = np.append(x0, 1.0+1.0/float(bits))
y = np.append(y0, 1.0+1.0/float(bits))

xg, yg = np.meshgrid(x[:],y[:])

h = np.loadtxt('h81.txt')
psi = np.loadtxt('psiMat81.txt')

i=5
#wut = u0[i*len(y):((i)*len(y)+len(x)),:]



#### PLOTTING ####
# cm.Spectral_r


fig=plt.figure()

count = 1
i = 3100
ax1=fig.add_subplot(1,2,1, aspect='equal')


#h = u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
#print u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:].shape
#psi = psi[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
#print psi[i*len(y)-i:((i)*len(y)+len(x))-i-1,:].shape
#ui = ut[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
#vi = vt[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]

#h = np.concatenate((h, np.ones((1,bits))))

#p = plt.pcolor(xg,yg,h,cmap=cm.Spectral_r,edgecolor='none')
#p.set_clim([np.min(u0[1*len(y)-1:((1)*len(y)+len(x))-1,:]), np.max(u0[1*len(y)-1:((1)*len(y)+len(x))-1,:])])


levels = np.arange(-20.0,21.0,1.0)


print "hm"
print y0.shape
print x0.shape
print psi.shape
print h.shape

CS = plt.contour(x0, y0, -psi, levels, colors='k', linewidths=np.array([.5]))
CS.set_clim([np.min(psi), np.max(psi)])
levels2 = [-20, -17, -14, -11, -8, -5]
plt.clabel(CS, levels2, fontsize=12,inline=1, fmt='%1.1f')
#plt.quiver(x0,y0,ui,vi)
plt.title("streamlines",fontsize=18)

ax1=fig.add_subplot(1,2,2, aspect='equal')
levels = np.arange(-0.5,0.5,.05)
#levels = [-.25, 0.0, .25]
print y0.shape
print x0.shape
print psi.shape
print h.shape
levels = np.arange(-.5,.5,.05)
IS = plt.contour(x0, y0, h, levels, colors='k', linewidths=np.array([.5]))
IS.set_clim([np.min(h), np.max(h)])
plt.clabel(IS, fontsize=12, manual='True',inline=1, fmt='%1.2f')
plt.yticks([],{'fontsize':'small'})

plt.title("isotherms",fontsize=18)
plt.xticks([])


count = count + 1
plt.xlim(x[0], x[-2])
plt.ylim(y[0], y[-2])

    
plt.subplots_adjust(bottom=.0, left=.05, right=.95, top=1.0, hspace=.3)

#cax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
#cbar = plt.colorbar(p, ticks = [int(np.min(u0[1*len(y)-1:((1)*len(y)+len(x))-1,:])),int(np.max(u0[1*len(y)-1:((1)*len(y)+len(x))-1,:]))],cax=cax,orientation='horizontal')


#cbar.set_label(r'TEMPERATURE',fontsize=10)

plt.savefig('r0718.png')


