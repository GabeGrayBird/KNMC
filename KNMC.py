import numpy as np
import random as r
import sys
"""
The goal of this simulation is to generate a random possible path for a photon from a given source based on the geometries of a given setup. Those angles will then be fed into the 
Klein-Nishina formula for a certain number of events to see what we should expect given our
setup.

There will be two main events that determine scattering angles. The intial pass through the aperture and the orientation of scattering in the target.
"""
h=4.135667e-15
me=.551
Ep=667.1e3
c=3e8
#eps=Ep/(511e3)
#print(eps)
re=2.82e-15
w1=300e-12
w2=600e-12
def gDg(ang, eps): #gamma' divided by gamma
    return (1/(1+(eps/511e3)*(1-np.cos(ang))))
def Tomp(ang):
    return .5*(re**2)*(1+(np.cos(ang))**2)
"""
def he(ang):
    return .5*(re**2)/(eps*(1-(np.cos(ang))))
"""
def KN(ang, eps):
    return .5*(re**2)*(gDg(ang,eps)**2)*(gDg(ang,eps)+(1./gDg(ang,eps))-(np.sin(ang))**2)

def randAng(minAng,maxAng): # return a random angle within the given range (in degrees)
    return r.uniform(minAng,maxAng)

def pos(g1, vlabg, t):
    return g1 + vlabg*t

def thetaKN(v1,v2): #dot product for getting the actual KN angle for count weight
    temp = np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    return temp
    
def quadSol(a, b, c): #Quadratic solver, used for finding time t which solves for hit position
    return np.array([(-b+np.sqrt(b**2-4*a*c))/(2*a), (-b-np.sqrt(b**2-4*a*c))/(2*a)])
#Coordinate system
#phi is the azimuthal angle (rotational/cone)
#theta is the beam angle
def scatter(thetaMax, L, R):
    thetaM=np.deg2rad(thetaMax)
    theta1=r.uniform(0, thetaM)
    phi1=r.uniform(0, 2*np.pi)
    #z is the beam line, and x and y are set up to make 
    #a right hand coordinate system
    x1=L*np.tan(theta1)*np.sin(phi1)
    y1=L*np.tan(theta1)*np.cos(phi1)
    z1=0
    #same process again but with different ranges
    thetaLab=r.uniform(0, np.pi)
    phiLab=r.uniform(0,2*np.pi)
    #lab velocity vector off the target 
    vLab=np.array([np.sin(thetaLab)*np.cos(phiLab), np.sin(thetaLab)*np.sin(phiLab), np.cos(thetaLab)])
    v1=np.array([x1, y1, L])
    a=vLab[0]**2+vLab[1]**2
    b=2*(x1*vLab[0] + y1*vLab[1])
    c=x1**2+y1**2 - R**2
    bothT=quadSol(a,b,c)
    t=bothT[0]
    if (t<=0): #if this solution is negative, use other one that is positive
        t=bothT[1]
    xFinal=x1+vLab[0]*t
    yFinal=y1+vLab[1]*t
    zFinal=z1+vLab[2]*t
    thetaFinal=thetaKN(vLab, v1) #the incident beam frame scattering angle (KN)
    temp = np.array([xFinal, yFinal, zFinal, thetaFinal])
    return temp

#define a range to  start and end at with steps, how many events
#rDet=radius of the detector, R=Radius of the chamber, L=distance from source
#to target, appAng=maximum angle determined from the appeture
def RunExp(degStart, degEnd, degStep, events, rDet, R, L,appAng):
    counts=[]
    bigSum=0
    degVals=np.arange(degStart, degEnd+degStep, degStep)
    for a in degVals:
        bigSum=0
        print(f'Doing {a} degree tests...',file = sys.stderr)
        for i in range(0, events, 1):
            if((i)%(events/10)==0):
                print(f'{100*(i/events)}%')
                
            data = scatter(np.deg2rad(appAng), L, R) #do a scattering event
            #data is structured [xFinal, yFinal, zFinal, thetaKN]
            #thetaHit=np.arctan(data[0]/data[2]) # arctan(x/z)
            thetaHit=np.arctan2(data[0],data[2]) # arctan(x/z)
            if((R*(thetaHit-np.deg2rad(a)))**2 + data[1]**2 < rDet**2):
                bigSum+=1*KN(data[3], 667.1e3) #weighing the count
                
        counts.append(bigSum) #append the sum for this range once done
    comb = np.column_stack((degVals, counts)) #format data
    return comb #return a 2D np array 

def RunExpOneAng(ang, events):
    rDet=.0171
    R=.508
    L=.1229
    appAng=5.08
    bigSum=0
    print(f'Doing {ang} degree tests...',file = sys.stderr)
    for i in range(0, events, 1):
        if((i)%(events/10)==0):
            print(f'{100*(i/events)}%',file = sys.stderr)
                
        data = scatter(np.deg2rad(appAng), L, R) #do a scattering event
        #data is structured [xFinal, yFinal, zFinal, thetaKN]
        #thetaHit=np.arctan(data[0]/data[2]) # arctan(x/z)
        thetaHit=np.arctan2(data[0],data[2]) # arctan(x/z)
        if((R*(thetaHit-np.deg2rad(ang)))**2 + data[1]**2 < rDet**2):
            bigSum+=1*KN(data[3], 667.1e3) #weighing the count
            
    #counts.append(bigSum) #append the sum for this range once done
    #comb = np.column_stack(ang, bigSum) #format data
    #return comb #return a 2D np array 
    return bigSum

ang=float(sys.argv[1])

finalsum=RunExpOneAng(ang, int(1e7)) #Run the sim real L .1229

print(ang,finalsum)
#print(f'{finalList[0]} {finalList[1]}')


"""
#Write the data to a txt file
knmcFile=open("MCdataBig.txt", "w")
np.savetxt("MCdataBig.txt", finalList)
knmcFile.close()
"""
