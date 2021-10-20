# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 14:52:56 2021

@author: Labor_NLO-Admin
"""
from NKTP_DLL import *
import Laser
import random
from time import sleep




def doRandomDisco():
    global laser
    i = 0
    while(i < 100):
        wavelength = random.randint(450,650)
        amplitude = random.randint(50,101)
        print(wavelength)
        laser.setLaserOutput(0x00, wavelength, amplitude)
        
        
        i = i+1
        sleep(0.01)
        
        
        
        
        
laser = Laser.Laser()
laser.openConnection()
laser.EmissionOn()
laser.resetInterLock()
sleep(1.0)


#laser.activateFSK()
laser.RFOn()
#laser.setLaserOutput(0x07, 510, 90)


for i in range(8):
    laser.setLaserOutput(i, 600, 0)


sleep(4.0)
laser.setLaserOutput(0x00, 650, 100)
print("Rot")
sleep(4.0)

doRandomDisco()
sleep(4.0)

for i in range(8):
    laser.setLaserOutput(i, 600, 0)

sleep(4.0)

#laser.deactivateFSK()
laser.closeConnection()
