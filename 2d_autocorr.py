# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:08:59 2022

@author: php20jo
"""

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.axes as ax
#import imageio as iio
from scipy.interpolate import interp1d
from scipy import interpolate



# def radial(matrix):
#     s=matrix.shape
#     dim=len(s)
#     center=[]
#     for i in range(len(s)):
#     	if  s[i] % 2 == 0:
#     		center.append(s[i]/2)
#     	else:
#     		center.append(s[i]/2+1)
#     if dim==2:
#     	ys, xs= np.meshgrid(np.arange(matrix.shape[1]),
#     	np.arange(matrix.shape[0]))
#     	rs=np.sqrt((ys-center[0])**2+(xs-center[1])**2)
#     # elif dim==3:
#     # 	ys, xs, zs = np.meshgrid(np.arange(matrix.shape[0]),
#     # 	np.arange(matrix.shape[1]),
#     # 	np.arange(matrix.shape[2]))
#     # 	rs=np.sqrt((ys-center[0])**2+(xs-center[1])**2+(zs-center[2])**2)
#     else:
#     	print('dimension error')
#     	return
    	
#     rmax=int(round(np.max(rs)))
#     # print(rs.shape)
#     # print(matrix.shape)
#     g, bins = np.histogram(rs, weights=matrix,bins=rmax)
#     counts, bins = np.histogram(rs,bins=rmax)
#     g /= counts
#     bins=(bins[1:]+bins[:-1])/2
#     #bins=bins*.025
#     return g, bins


# def projections_3d(matrix):
# 	#make the vectors for each G function 
# 	ycor=np.arange(0,ly,1)
# 	xcor=np.arange(0,lx,1)
# 	zcor=np.arange(0,lz,1)
# 	#Calculate the corelation length Eq.(1.6) Anderson book
# 	cly=2*np.sum(matrix[:,lx/2,lz/2])
# 	clx=2*np.sum(matrix[ly/2,:,lz/2])
# 	clz=2*np.sum(matrix[ly/2,lx/2,:])
# 	#Each projection y,x,z based on Eq.(1.7) Anderson book
# 	Gy=[]
# 	for z in range(len(zcor)):	
# 		Gy.append((2/cly)*np.sum(matrix[:,lx/2,zcor[z]]))
# 	Gx=[]
# 	for y in range(len(ycor)):
# 		Gx.append((2/clx)*np.sum(matrix[ycor[y],:,lz/2]))
# 	Gz=[]
# 	for x in range(len(xcor)):	
# 		Gz.append((2/clz)*np.sum(matrix[ly/2,xcor[x],:]))
# 	return ycor, Gy, xcor, Gx, zcor, Gz, cly, clx, clz 

def projections_2d(matrix):
    #make the vectors for each G function 
    ycor=np.arange(0,ly,1)
    xcor=np.arange(0,lx,1)
    # print(lx/2)
    #Calculate the corelation length Eq.(1.6) Anderson book
    cly=2*np.sum(matrix[:,int(lx/2)])
    clx=2*np.sum(matrix[int(ly/2),:])
    #Each projection y,x,z based on Eq.(1.7) Anderson book
    Gy_x=[]
    for x in range(len(xcor)):	
    	Gy_x.append((2/cly)*np.sum(matrix[:,xcor[x]]))
    Gx_y=[]
    for y in range(len(ycor)):
    	Gx_y.append((2/clx)*np.sum(matrix[ycor[y],:]))
    
    return ycor, Gy_x, xcor, Gx_y, cly, clx

# def proj_r(gamma):
#     gamma_r,radii=radial(gamma)
#     plt.figure()
#     plt.plot(radii,gamma_r)
#     gamma_r=[a/np.max(gamma_r) for a in gamma_r]
#     spline = interpolate.splrep(radii, gamma_r, s=0)
#     delta_r=0.5
#     rnew=np.arange(np.min(radii),np.max(radii),delta_r)
#     zcor=rnew
#     plt.plot(zcor,interpolate.splint(zcor[0],zcor[-1],spline))
#     clr=2*interpolate.splint(zcor[0],zcor[-1],spline)
#     print(clr)	
#     Gr=[]
#     for z in range(len(zcor)):
#     	term1=interpolate.splev(zcor[-1], spline, der=0)*np.sqrt(zcor[-1]**2-zcor[z]**2)
#     	zval=zcor[z:]
#     	term2=delta_r*np.sum(np.sqrt(zval**2-zval[0]**2)*interpolate.splev(zval, spline, der=1))
#     	Gr.append(2/clr*(term1-term2))
    # return Gr,zcor


def autocorr(matrix):
	# Calculate the autocorrelation function of tomography section
	f = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(matrix)))	
	P = np.abs(f)**2
	transform = abs(np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(P))))
	return transform

#%%

scale=64 #pix/um

test = plt.imread('sio2_002.002_th.png')
#test=iio.imread('sio2_002.002_th.png')
test=test[:,:-1]
print(test.shape)

plt.figure() 
plt.imshow(test,cmap='gray')
plt.title('raw image')
plt.gca().set_xscale('linear')

#%%


s=test.shape  
ly=s[0]
lx=s[1]
# lz=s[2]
dely=1;
delx=1;
# delz=1;
N=lx*ly#*lz
N0=np.sum(test)
phi_n=N0/N
phi_v=phi_n
V=N*delx*dely#*delz

gamma=autocorr(test)

plt.figure()
plt.imshow(gamma)
plt.title('autocorr')

ycor, Gy_x, xcor, Gx_y, cly, clx=projections_2d(gamma)


plt.figure()
plt.plot(xcor[:int(len(xcor)/2)],Gy_x[int(len(xcor)/2):-1],'r',label='Gy(x)')
plt.legend()
plt.figure()
plt.plot(ycor[:int(len(ycor)/2)],Gx_y[int(len(ycor)/2):],'b--',label='Gx(y)')
plt.legend()












