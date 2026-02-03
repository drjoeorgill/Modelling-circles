# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:09:10 2023

@author: php20jo
"""

import numpy as np
import matplotlib.pyplot as plt
#from random import randint
#from scipy.spatial import distance 
#from mpl_toolkits.mplot3d import Axes3D
#from itertools import product
import imageio as iio
import os


# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# Change the current working directory
os.chdir('C:/Users/pczjo/OneDrive - The University of Nottingham/Desktop/Modeling_circles')

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))


def autocorr(matrix):
	# Calculate the autocorrelation function of tomography section
	f = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(matrix)))	
	P = np.abs(f)**2
	transform = abs(np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(P))))
	return transform

def Fourier(matrix):
	# Calculate the fourier transfom of an image 
    f = abs(np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(matrix))))
    return f

def Fourier2(matrix):
	# Calculate the fourier transform squared of an image
    f = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(matrix)))	
    P = np.abs(f)**2
    return P

def inFourier(matrix):
	# Calculate the 	inverse fourier of an image
    transform = abs(np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(matrix))))
    return transform

# i havent made a function yet that does a 2d projection of the fourier

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

#%% raw image

#test=iio.imread('test_circles.png')
#test=test[:,:160]
#print(test.shape)

test = plt.imread('test_circles.png')
test=test[:,:160]
print(test.shape)


plt.figure() 
plt.imshow(test,cmap='gray')
plt.title('raw image')
plt.gca().set_xscale('linear')

#%% Fourier

f = Fourier(test)
plt.figure()
plt.imshow(np.log(f))
plt.colorbar()
plt.title('Fourier')

#%% autocorr

gamma = autocorr(test)
plt.figure()
plt.imshow(gamma)
plt.title('autocorr')



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

F=Fourier(test)

plt.figure()
plt.imshow(np.log(F))
plt.title('Fourier')
#%%
ycor, Gy_x, xcor, Gx_y, cly, clx=projections_2d(F)


plt.figure()
plt.semilogy(xcor[:int(len(xcor)/2)],Gy_x[int(len(xcor)/2):],'r',label='Gy(x)')
plt.legend()
plt.figure()
plt.plot(ycor[:int(len(ycor)/2)],Gx_y[int(len(ycor)/2):],'b--',label='Gx(y)')
plt.legend()

#%% 2d fourier

'''
#this needs to use a new function that intergrates along a line from the center
yf, Gy_x, xf, Gx_y, cly, clx=projections_2d(f)


plt.figure()
plt.semilogy(xcor[:int(len(xcor)/2)],Gy_x[int(len(xcor)/2):],'r',label='Gy(x)')
plt.legend()
plt.figure()
plt.plot(ycor[:int(len(ycor)/2)],Gx_y[int(len(ycor)/2):],'b--',label='Gx(y)')
plt.legend()
'''