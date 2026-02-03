# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 15:11:49 2022

@author: php20jo
"""

import numpy as np
import matplotlib.pyplot as plt
from random import randint
from scipy.spatial import distance 
from mpl_toolkits.mplot3d import Axes3D
from itertools import product

def autocorr(matrix):
	# Calculate the autocorrelation function of tomography section
	f = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(matrix)))	
	P = np.abs(f)**2
	transform = abs(np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(P))))
	return transform

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

def ordered_spheres(sides, radius, x_separation, y_separation):
    '''
    Parameters
    ----------
    sides : int
        size of matrix.
    radius : int
        radius of each circle.
    x_separation : int
        separation in x axis.
    y_separation : int
        separation in y axis.

    Returns
    -------
    a matrix of circles in a square packing pattern with.

    '''
    box = np.ones((sides,sides))
    #define x coordinates for circle centers
    xcenters = [radius] #first circle can be pushed up against the side of the box
    i=1
    while max(xcenters)<(sides-radius):
        xcenters.append(radius+x_separation*i)
        i+=1
    if max(xcenters)>(sides-radius):
        xcenters.pop() #remove any circles that will overlap sides
    #define y coordinates for circle centers
    ycenters = [radius]
    j=1
    while max(ycenters)<(sides-radius):
        ycenters.append(radius+y_separation*j)
        j+=1
    if max(ycenters)>(sides-radius):
        ycenters.pop()
    
    #define all points in a circle around each coordinate as 0
    for m in xcenters:
        for n in ycenters:
            bounds = [m-radius, m+radius, n-radius, n+radius]    
            for i in range(bounds[0],bounds[1]):
                for j in range(bounds[2],bounds[3]):
                        xt = i - m
                        yt = j - n
                        if yt**2+xt**2<radius**2:
                            box[i,j]=0
    return(box)

def hcp_circles(sides, radius, x_separation, y_separation):
    '''
    Parameters all in nm
    ----------
    sides : int
        size of matrix.
    radius : int
        radius of each circle.
    x_separation : int
        separation in x axis divided by 2.
    y_separation : int
        separation in y axis divided by 2.

    Returns
    -------
    a matrix of circles in a square packing pattern with.

    '''
    box = np.ones((sides,sides))
    #define x coordinates for circle centers
    xcenters1 = [radius] #first circle can be pushed up against the side of the box
    i=1
    while max(xcenters1)<(sides-radius):
        xcenters1.append(radius+x_separation*2*i)
        i+=1
    if max(xcenters1)>(sides-radius):
        xcenters1.pop() #remove any circles that will overlap sides
    #second lot of x coordinates
    xcenters2 = [radius+x_separation] #first circle can be pushed up against the side of the box
    i=1
    while max(xcenters2)<(sides-radius):
        xcenters2.append(radius+x_separation+x_separation*2*i)
        i+=1
    if max(xcenters2)>(sides-radius):
        xcenters2.pop() #remove any circles that will overlap sides
    #define y coordinates for circle centers
    ycenters1 = [radius]
    j=1
    while max(ycenters1)<(sides-radius):
        ycenters1.append(radius+y_separation*2*j)
        j+=1
    if max(ycenters1)>(sides-radius):
        ycenters1.pop()
    #second lot of y coordinates
    ycenters2 = [radius+ y_separation]
    j=1
    while max(ycenters2)<(sides-radius):
        ycenters2.append(radius+y_separation+y_separation*2*j)
        j+=1
    if max(ycenters2)>(sides-radius):
        ycenters2.pop()
    
    
    #define all points in a circle around each coordinate as 0
    for m in xcenters1:
        for n in ycenters1:
            bounds = [m-radius, m+radius, n-radius, n+radius]    
            for i in range(bounds[0],bounds[1]):
                for j in range(bounds[2],bounds[3]):
                        xt = i - m
                        yt = j - n
                        if yt**2+xt**2<radius**2:
                            box[i,j]=0
    for m in xcenters2:
        for n in ycenters2:
            bounds = [m-radius, m+radius, n-radius, n+radius]    
            for i in range(bounds[0],bounds[1]):
                for j in range(bounds[2],bounds[3]):
                        xt = i - m
                        yt = j - n
                        if yt**2+xt**2<radius**2:
                            box[i,j]=0
    return(box)

#%%

co = ordered_spheres(100,5,20,15)


plt.figure()
plt.imshow(co,cmap='gray')

#%%

hcp = hcp_circles(1000,50,100,87)
plt.figure()
plt.imshow(hcp,cmap='gray')
plt.xlabel('10 nm')
plt.ylabel('10 nm')

#%%

s=hcp.shape  
ly=s[0]
lx=s[1]
# lz=s[2]
dely=1;
delx=1;
# delz=1;
N=lx*ly#*lz
N0=np.sum(hcp)
phi_n=N0/N
phi_v=phi_n
V=N*delx*dely#*delz

gamma=autocorr(hcp)
#%%
plt.figure()
plt.imshow(gamma)
#%%
ycor, Gy_x, xcor, Gx_y, cly, clx=projections_2d(gamma)


plt.figure()
plt.plot(xcor[:int(len(xcor)/2)],Gy_x[int(len(xcor)/2):],'r',label='Gy(x)')
plt.legend()
plt.figure()
plt.plot(ycor[:int(len(ycor)/2)],Gx_y[int(len(ycor)/2):],'b',label='Gx(y)')
plt.legend()
    
