# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 11:05:27 2021

@author: Labor_NLO-Admin
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import Pattern as pt

def getSineIntegerLayout(startWave, endWave, step, amplitude, yOffset, phase, periode):
    frequ = 2*np.pi/periode
    area = list(range(startWave, endWave+step, step))
    ar = np.array(area)
    function = amplitude*np.sin(ar*frequ+phase)+yOffset
    #function = np.floor(function)
    #function *= 1000
    #function = function.astype(int)
    return (ar, function)
                         

pattern = getSineIntegerLayout(470, 660, 1, 0.5, 0.5, 0, 184*2)

plt.plot(pattern[0],pattern[1])
plt.show()

periods = np.linspace(184*2,3,184*4)
print(periods)




for i in range(184):
            #set laser output
            self.laser.setLaserOutputPM(0, int(pattern.wavelengths[i]*1000), int(pattern.intensities[i]/corr_func(pattern.wavelengths[i])))
            #append wavelength
            plotWave.append(pattern.wavelengths[i])
            #append intensity
            plotInt.append(int(pattern.intensities[i]))
            #get wavelengths of spectrometer
            current_wave = self.spec.wavelengths()
            #get current spectrum
            current_int = self.spec.intensities()
            #get image
            img = self.cam.takeImage()
            #save spectrum
            self.stsToCSV(current_int,current_wave,r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\pattern_test\\"+str(pattern.wavelengths[i])+"nm_STS.spec")
            #save image
            self.imageToCSV(img, r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\pattern_test\\"+str(pattern.wavelengths[i])+"nm_CEMOS.data")