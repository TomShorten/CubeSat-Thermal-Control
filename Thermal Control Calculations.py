# importing librarys
import numpy as np
import math as m
import matplotlib
from matplotlib import pyplot as plt

# Defining our orbit constants

a = 1150000                     # Orbit altitude
i = 86                          # Inclination
e = 0                           # eccentricity
mu = 4.90486959*m.pow(10,12)    # Moon gravitational constant
rmoon = 1737100                 # Moon radius
T_max = 313.15                  # Upper T limit
T_min = 283                     # Lower T limit
T_space = 2.7                   # Space temp
sigma = 5.67*10**(-8)           # Stefan-Boltzmann constant

# More calculated constants 

rorb = rmoon + a
period = 2*m.pi*m.sqrt(m.pow(rorb,3)/mu)
radConv = m.pi/180
degConv = 180/m.pi 
eclipseAngleLower = m.pi - ( m.pi/2 - m.acos(rmoon/rorb))
eclipseAngleUpper = m.pi + ( m.pi/2 - m.acos(rmoon/rorb))
eclipseTime = period * ((eclipseAngleUpper-eclipseAngleLower) / (2*m.pi))

# Display the orbit time (s) and eclipse time (s) in console

print(round(period), 'seconds per orbit')
print(round(eclipseTime), 'seconds of eclipse per orbit')
print("")

# Define CubeSat face area

Face_1 = ( 0.1 * 0.2 ) # SP alpha = 0.67 unchangeable 
Face_2 = ( 0.2 * 0.3 ) # MLI alpha depends on No. layers
Face_3 = ( 0.2 * 0.3 ) # MLI alpha depends on No. layers
Face_4 = ( 0.1 * 0.2 ) # Solar absorber, high alpha, low e 0.3/0.02
Face_5 = ( 0.1 * 0.3 ) # Biomechanical louver alpha = 0.57, e = 0.03
Face_6 = ( 0.1 * 0.3 ) # Biomechanical louver

# Define CubeSat face optical properties

alphaF1 = 0.25      # Gimble / Solar Array
alphaMLI = 0.14     # 2x3 face
alphaF4 = 0.25      # Paylaod face
alphaF56 = 0.25     # Side Faces
eF1 = 0.33          # Gimle / Solar Array
eMLI = 0.035        # 2x3 face
eF4 = 0.33          # Payload Face
eF56 = 0.33         # Side Faces

# Define constants for heat balance equations

sin4 = m.sin(radConv*4)                     # Inclination angle effect
cos4 = m.cos(radConv*4)                     # Inclination angle effect
theta = np.linspace(0,360,360,False)        # Create an array of true anomalies
beta = m.cos(radConv*4)                     # beta angle for moon heat input
Jsun = 1361                                 # Sun heat input
C1 = 1335                                   # C1 from moon heat equation
C2 = 5                                      # C2 from moon heat equation

# Create array for sun and moon heat input

Q_inSun = []
Q_inMoon = []

# Itterative Sun heat input for each degree in orbit


for i in theta:
    if 0 <= i <= (degConv*eclipseAngleLower):
        if 0 <= i <= 90:
            Q1 = Face_1*m.cos(radConv*i)*alphaF1*Jsun
            Q4 = 0
        else:
            Q1 = 0
            Q4 = -Face_4*m.sin(radConv*i + (m.pi/2))*alphaF4*Jsun
        Q2 = Face_2*m.sin(radConv*i)*alphaMLI*Jsun
        Q3 = 0
    
    elif  degConv*eclipseAngleLower < i <= degConv*eclipseAngleUpper:
        Q1 = Q2 = Q3 = Q4 = 0
    
    elif degConv*eclipseAngleUpper < i <= 360:
        if degConv*eclipseAngleUpper < i <= 270:
            Q1 = 0
            Q4 = -Face_4*m.sin(radConv*i + (m.pi/2))*alphaF4*Jsun
        else:
            Q1 = Face_1*m.cos(radConv*i)*alphaF1*Jsun
            Q4 = 0
        Q2 = 0
        Q3 = Face_3*m.cos(radConv*i + (m.pi/2))*Jsun*alphaMLI
        
    Q_inSun.append(Q1+Q2+Q3+Q4)

# Itterative Moon heat input for each degree in orbit

for i in theta:
    if 0 <= i <= 90:
        Q4 = ((C1 - C2)*beta*m.cos(radConv*i)+C2)*Face_4*alphaF4
        
    elif 90 < i <= 270:
        Q4 = C2*Face_4*alphaF4
        
    elif 270 < i <= 360:
        Q4 = ((C1 - C2)*beta*m.cos(radConv*i)+C2)*Face_4*alphaF4
    
        
    Q_inMoon.append(Q4)


# Sum of Moon and Sun heat input

Q_total = [Q_inSun[i]+ Q_inMoon[i] for i in range(len(Q_inSun))]

# Initialising temperature array, heater power, heater input

Temp = []
Heater = np.linspace(0,20,100,False)
HeaterInput = []
Tempnoheater = []
Tlimit = []
AeSum = (Face_1*eF1+Face_2*eMLI+Face_3*eMLI+Face_4*eF4+Face_5*eF56+Face_6*eF56)
SigAe = sigma*AeSum
# Itterative calculation of Temperature through 0 deg ascension #### orbit 

for i in theta:                                     
    if 0 <= i < 90:
        Q1 = Face_1*m.cos(radConv*i)*alphaF1*Jsun
        Q4 = 0
        Q4MOON = ((C1 - C2)*beta*m.cos(radConv*i)+C2)*Face_4*alphaF4
        Q2 = Face_2*m.sin(radConv*i)*alphaMLI*Jsun
        Q3 = 0
        Tto4 = T_space**4 + ( Q1 + Q2 + Q3 + Q4 + Q4MOON )/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        z = 0
        for j in Heater:
            if T >= T_min:
                break
            else:
                Tto4 = T_space**4 + ( j + Q1 + Q2 + Q3 + Q4 + Q4MOON )/(SigAe)
                T = Tto4**(1/4)
                z = j

                   
       

    elif 90 <= i <= (degConv*eclipseAngleLower):
        Q1 = 0
        Q4 = -Face_4*m.sin(radConv*i + (m.pi/2))*alphaF4*Jsun*cos4
        Q4MOON = C2*Face_4*alphaF4
        Q2 = Face_2*m.sin(radConv*i)*alphaMLI*Jsun
        Q3 = 0
        Tto4 = T_space**4 + (Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        z = 0
        for j in Heater:
            if T >= T_min:
                break
            else:
                Tto4 = T_space**4 + ( j + Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
                T = Tto4**(1/4)
                z = j
        
        
        
    elif  degConv*eclipseAngleLower < i <= degConv*eclipseAngleUpper:
        Q1 = Q2 = Q3 = Q4 = 0
        Q4MOON = C2*Face_4*alphaF4
        Tto4 = T_space**4 + (Q1 + Q2 + Q3 + Q4 + Q4MOON )/(SigAe)
        T = Tto4**(1/4)
        Tto4nh = T_space**4 + ( Q1 + Q2 + Q3 + Q4 + Q4MOON )/(SigAe)
        Tnoheater = Tto4nh**(1/4)
        z = 0
        for j in Heater:
            if T >= T_min:
                break
            else:
                Tto4 = T_space**4 + ( j + Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
                T = Tto4**(1/4)
                z = j
    
    elif degConv*eclipseAngleUpper < i <= 270: 
        Q1 = 0
        Q4 = -Face_4*m.sin(radConv*i + (m.pi/2))*alphaF4*Jsun*cos4
        Q4MOON = C2*Face_4*alphaF4
        Q2 = 0
        Q3 = Face_3*m.cos(radConv*i + (m.pi/2))*Jsun*alphaMLI
        Tto4 = T_space**4 + (Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        z = 0
        for j in Heater:
            if T >= T_min:
                break
            else:
                Tto4 = T_space**4 + ( j + Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
                T = Tto4**(1/4)
                z = j
        
    elif 270 < i <= 360:
        Q1 = Face_1*m.cos(radConv*i)*alphaF1*Jsun*cos4
        Q4 = 0
        Q4MOON = ((C1 - C2)*beta*m.cos(radConv*i)+C2)*Face_4*alphaF4
        Q2 = 0
        Q3 = Face_3*m.cos(radConv*i + (m.pi/2))*Jsun*alphaMLI
        Tto4 = T_space**4 + (Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        z = 0
        for j in Heater:
            if T >= T_min:
                break
            else:
                Tto4 = T_space**4 + ( j + Q1 + Q2 + Q3 + Q4 + Q4MOON)/(SigAe)
                T = Tto4**(1/4)
                z = j
        
    Temp.append(T)
    HeaterInput.append(z)
    Tempnoheater.append(Tnoheater)
    Tlimit.append(T_max)

# Displaying maximum and minimum temperature through orbit, and max heater input

print(round(np.min(Temp), 2))
print(round(np.max(Temp), 2))
print(round(np.max(HeaterInput), 2))

# Initilise new array for new orbit

Q_inSun_90 = []
Q_inMoon_90 = []

# Itterative Sun heat input for each degree in 90 ascension #### orbit

for i in theta:
    if 0 <= i <= 180:
        if 0 <= i <= 90:
            Q2 = Face_2*m.cos(radConv*i)*alphaMLI*Jsun*sin4
            Q3 = 0
        else:
            Q2 = 0
            Q3 = Face_3*m.sin(radConv*i - (m.pi/2))*Jsun*alphaMLI*sin4
        Q4 = 0
        Q1 = Face_1*m.sin(radConv*i)*alphaF1*Jsun*sin4
        Q5 = Face_5*alphaF56*Jsun*cos4

    elif 180 < i <= 360:
        if 180 < i <= 270:
            Q2 = 0
            Q3 = Face_3*m.sin(radConv*i - (m.pi/2))*Jsun*alphaMLI*sin4
        else:
            Q2 = Face_2*m.sin(radConv*i - (3/2)*m.pi)*alphaMLI*Jsun*sin4
            Q3 = 0
        Q1 = 0
        Q4 = Face_4*m.sin(radConv*i - m.pi)*alphaF4*Jsun*sin4
        Q5 = Face_5*alphaF56*Jsun*cos4
    Q_inSun_90.append(Q1+Q2+Q3+Q4+Q5)

# Itterative Moon heat input for each degree in 90 ascension #### orbit

for i in theta:
    Q4 = C2*Face_4*alphaF4

    Q_inMoon_90.append(Q4)

# Total heat input calcualation

Q_total_90 = [Q_inSun_90[i]+ Q_inMoon_90[i] for i in range(len(Q_inSun_90))]

# Initialise new temperature during orbit array

Temp_90 = []
HeaterInput_90 = []
Tempnoheater_90 = []

# Itterative calculation of Temperature through 90 degree ascention #### orbit

for i in theta:                                     
    if 0 <= i < 180:
        if 0 <= i <= 90:
            Q2 = Face_2*m.cos(radConv*i)*alphaMLI*Jsun*sin4
            Q3 = 0
        else:
            Q2 = 0
            Q3 = Face_3*m.sin(radConv*i - (m.pi/2))*Jsun*alphaMLI*sin4
        Q4 = 0
        Q1 = Face_1*m.sin(radConv*i)*alphaF1*Jsun*sin4
        Q5 = Face_5*alphaF56*Jsun*cos4
        Q4MOON = C2*Face_4*alphaF4
        Tto4 = T_space**4 + ( Q1 + Q2 + Q3 + Q4 + Q5 + Q4MOON)/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        
                   
       

    elif 180 <= i <= 360:
        if 180 <= i <= 270:
            Q2 = 0
            Q3 = Face_3*m.sin(radConv*i - (m.pi/2))*Jsun*alphaMLI*sin4
        else:
            Q2 = Face_2*m.sin(radConv*i - (3/2)*m.pi)*alphaMLI*Jsun*sin4
            Q3 = 0
        Q1 = 0
        Q4 = Face_4*m.sin(radConv*i - m.pi)*alphaF4*Jsun*sin4
        Q5 = Face_5*alphaF56*Jsun*cos4
        Q4MOON = C2*Face_4*alphaF4
        Tto4 = T_space**4 + (Q1 + Q2 + Q3 + Q4 + Q4MOON + Q5)/(SigAe)
        T = Tto4**(1/4)
        Tnoheater = Tto4**(1/4)
        
    Temp_90.append(T)
   # HeaterInput_90.append(z)
    Tempnoheater_90.append(Tnoheater)


# Showing all graphs


# Heat input. total, sun and moon

matplotlib.rcParams.update({'font.size': 18})
plt.figure(1, figsize=(8,8))
plt.plot(theta, Q_inSun, color='lavenderblush')
plt.plot(theta, Q_inMoon, color='palegreen')
plt.plot(theta, Q_total, color='skyblue')
plt.xlabel('True Anomaly (-180)')
plt.ylabel('Total heat from all faces (W)')
plt.title('Total heat input')
plt.show()

# 0 deg ascention #### temp during orbit

matplotlib.rcParams.update({'font.size': 18})
plt.figure(1, figsize=(8,8))
plt.plot(theta, Temp, color='blue')
plt.plot(theta, Tempnoheater, color='red')
plt.plot(theta, Tlimit, color='red')
plt.xlabel('True Anomaly (-180)')
plt.ylabel('Temperature K')
plt.title('Temperature variation')
plt.show()

# Heater power required during orbit

matplotlib.rcParams.update({'font.size': 18})
plt.figure(1, figsize=(8,8))
plt.plot(theta, HeaterInput, color='blue')
plt.xlabel('True Anomaly (-180)')
plt.ylabel('Heat input W')
plt.title('Heater power during orbit')
plt.show()

# Sun, moon and total heat input 90 deg ascention

matplotlib.rcParams.update({'font.size': 18})
plt.figure(1, figsize=(8,8))
plt.plot(theta, Q_inSun_90, color='skyblue')
plt.plot(theta, Q_total_90, color='green')
plt.xlabel('True Anomaly (-180)')
plt.ylabel('Total heat in (W)')
plt.title('Sun heat input')
plt.show()

# 90 deg ascention #### temp during orbit

matplotlib.rcParams.update({'font.size': 18})
plt.figure(1, figsize=(8,8))
plt.plot(theta, Temp_90, color='blue')
plt.plot(theta, Tempnoheater_90, color='red')
#plt.plot(theta, Tlimit, color='red')
plt.xlabel('True Anomaly (-180)')
plt.ylabel('Temperature K')
plt.title('Temperature variation')
plt.show()