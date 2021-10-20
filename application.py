#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 18:09:21 2021

@author: Jonas Brandhoff
"""
import scipy.interpolate as itp
import scipy.optimize as opt
import os
import sys
from NKTP_DLL import *
import Laser
import random
from PyQt5.QtWidgets import (
    QApplication, QDialog, QMainWindow, QMessageBox
)
from PyQt5.uic import loadUi
from PyQt5.QtWidgets import QDialog, QFileDialog
from datetime import datetime
from MainWindow import Ui_MainWindow
import pandas as pd
import numpy as np
import time
import asyncio
from PyQt5.QtCore import QTimer
from asyncqt import QEventLoop, asyncSlot
import Camera.camera as ca
import IPython
import Pattern as pt
from seabreeze.spectrometers import Spectrometer as sm
import matplotlib.pyplot as plt
import time
#IPython.get_ipython().run_line_magic('matplotlib', 'qt')


class Window(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.connectSignalsSlots()
        self.EmissionOn = False
        self.RFOn = False
        self.laser = Laser.Laser()
        self.laser.openConnection()
        
        self.timer=QTimer()
        self.runningTimer = False
        self.cam = ca.Camera()
        self.spec = sm.from_first_available()
        self.PatternHelper = pt.PatternHelper()
        
        
        dataCalibration = pd.read_csv('Calibrations/LaserSTSCalibration_merged.spec', sep=',')
        self.caliSTSint = dataCalibration['intensity'].values
        self.caliSTSwave = dataCalibration['wavelength'].values
        
        dataCalibration = pd.read_csv('Calibrations/LaserCMOSCalibration.spec', sep=',')
        self.caliCMOSint = dataCalibration['intensity'].values
        self.caliCMOSwave = dataCalibration['wavelength'].values

    """
    HIER WERDEN DIE EINZELNEN SIGNALS MIT DEN SLOTS VERBUNDEN
    siehe QT5 Doc
    """
    def connectSignalsSlots(self):
       self.btnExit.clicked.connect(self.exitOut)
       
       """
    
    LASER zeug
    
    """
       self.btnToggleEmission.clicked.connect(self.toggleEmission)
       self.btnToggleRF.clicked.connect(self.toggleRF)
       self.btnSetOutput.clicked.connect(self.changeLaser)
       self.btnRandomPattern.clicked.connect(self.freqSweep)
       """
       Power curve
       
       """
       self.btnTakeCurve.clicked.connect(self.Take2DCurve)
       
       
       """
       Camera stuff
       """
       self.pushButton.clicked.connect(self.camerShot)
       self.btnFitGaussian.clicked.connect(self.fitCurrentShot)
       self.pushButton_2.clicked.connect(self.roiCurrentShot)
    def toggleEmission(self):
        if self.EmissionOn:
            self.laser.EmissionOff()
            self.EmissionOn = False
            self.btnToggleEmission.setStyleSheet("border: 3px solid rgba(245,0,0,255);\n""border-radius: 10px;\n""background-color: rgba(245,0,0,255);")
            
            
        else:
            self.laser.resetInterLock()
            self.laser.EmissionOn()
            self.EmissionOn = True
            self.btnToggleEmission.setStyleSheet("border: 3px solid rgba(51, 181, 25,255);\n""border-radius: 10px;\n""background-color: rgba(51, 181, 25,255);")

    def toggleRF(self):
        if self.RFOn:
            self.laser.RFOff()
            self.RFOn = False
            self.btnToggleRF.setStyleSheet("border: 3px solid rgba(245,0,0,255);\n""border-radius: 10px;\n""background-color: rgba(245,0,0,255);")
            
            
        else:
            self.laser.RFOn()
            self.RFOn = True
            self.btnToggleRF.setStyleSheet("border: 3px solid rgba(51, 181, 25,255);\n""border-radius: 10px;\n""background-color: rgba(51, 181, 25,255);")

    

    def changeLaser(self):
        wav = self.spinWave.value()
        num = self.spinLaserNum.value()
        amp = self.spinAmp.value()
        self.laser.setLaserOutputPM(num, wav, amp)


 
        
    def exitOut(self):
       
        self.laser.RFOff()
        self.laser.EmissionOff()
        self.laser.closeConnection()
        self.cam.dispose()
        self.spec.close()
        os._exit(0)
        self.close()
        loop.stop()
        
    
    
    """
    --------------------------------------------------------------------------------------------------------
    --------------------------------GHOST PATTERNS----------------------------------------------------------
    --------------------------------------------------------------------------------------------------------
    """
    def startRandomPattern(self):
        
        if self.runningTimer:
            self.timer.stop()
            self.runningTimer = False
            self.timer.timeout.disconnect()
            return
        
        print("RANDOM PATTERN STARTED")
        self.timer.timeout.connect(self.doRandomPattern)
        self.runningTimer = True
        self.timer.start(1)
    def doRandomPattern(self):
        for i in range(8):
            wavelength = random.randint(4300,7000)
            amplitude = random.randint(1,101)
            self.laser.setLaserOutputPM(i, wavelength*100, amplitude*10)
    
    
    
    
    def imageToCSV(self, image, filename):
        df = pd.DataFrame()
        for i, im in enumerate(image):
            df.insert(0,"image"+str(i),im)
        df.to_csv(filename)
    
    def stsToCSV(self, intensity, wave, filename):
        df = pd.DataFrame()
        df.insert(0,"wavelength",wave)
        df.insert(0,"intensity",intensity)
        df.to_csv(filename)
    
    def getIndexForWave(self, wavelenght):
        exists = wavelenght in self.caliSTSwave
        if exists:
            index = list(self.caliSTSwave).index(wavelenght)
            return index
        else:
            print("ERROR WAVELENGTH NOT PRESENT RETURNING 0")
            return 0
    
    def getCalibratedIntensityFactor(self, wavelenght, STS):
        extraFunction = 1.5
        if STS:
            return 1/(extraFunction*self.caliSTSint[self.getIndexForWave(wavelenght)])
        else:
            return 1/self.caliCMOSint[self.getIndexForWave(wavelenght)]
    
    
    @asyncSlot()
    async def fastSweepThrouPattern(self):
        pattern = self.PatternHelper.getSineIntegerLayout(430,700,1, 500, 900, 0,0* 1/(730-430)*8*np.pi)
        plotWave = []
        plotInt = []
        self.spec.integration_time_micros(2000)
        for i in range(len(pattern.intensities)):
            self.laser.setLaserOutputPM(0, pattern.wavelengths[i]*1000, int(pattern.intensities[i]*self.getCalibratedIntensityFactor(pattern.wavelengths[i],True)))
            plotWave.append(pattern.wavelengths[i])
            plotInt.append(int(pattern.intensities[i]*self.getCalibratedIntensityFactor(pattern.wavelengths[i],True)))
            current_wave = self.spec.wavelengths()
            current_int = self.spec.intensities()
            self.stsToCSV(current_int,current_wave,"sinePattern/"+str(pattern.wavelengths[i])+"nm_STS.spec")
            schritteGes = len(pattern.wavelengths)
            momSchritt = i+1
            percent = momSchritt/schritteGes
            print(percent)
            self.progressCurve.setValue(percent*100)
        fig, ax = plt.subplots(1, gridspec_kw={'wspace': 0}, facecolor='w',figsize=(11,8))
        ax.plot(plotWave,plotInt, label='Spectrum')
        
        
    #function to run predefined patterns, still "hard coded"
    @asyncSlot()
    async def freqSweep(self):
        #define maximum frequency for 1 nm steps
        max_frequency = (450e-9)**2/(2*np.pi*2e-9)
        #define frequency array
        frequencies = np.linspace(0,max_frequency,184)
        #set integration time of spectrometer
        self.spec.integration_time_micros(5000)
        #set integration time of camera
        self.cam.setExposure(60)
        #set time out of camera
        self.cam.setTimeOut(300)
        
        pattern_list = ["frequency_sweep", "sine","cosine", "random"]
        #calibration of input intensity
        dataMatrix = pd.read_csv(r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\calibration\calibration_matrix_23_09.csv", sep=',')
        dataMatrix.drop(dataMatrix.columns[0], axis=1, inplace=True)
        np.set_printoptions(threshold=sys.maxsize)
        D = dataMatrix.to_numpy()
        correction = D[:,-1]
        correction /= np.min(correction)
        #define correction function of input intensities
        corr_func = itp.interp1d(np.linspace(470,660,39),correction)
        #set start time
        t0 = time.time()
        patterns = ["sine","cosine","sini","cosini"]
        #loop over all frequencies
        for k in range(len(frequencies)):
            print(k)
            if k != 20:
                continue
            #loop over sine, inverse sine, cosine and inverse cosine
            for l in patterns:
                if l == "sine":
                    pattern = self.PatternHelper.getSineIntegerLayout(470e-9, 654e-9, 1e-9, 0.5, 0.5, 0, frequencies[k])
                if l == "cosine":
                    pattern = self.PatternHelper.getSineIntegerLayout(470e-9, 654e-9, 1e-9, 0.5, 0.5, np.pi/2, frequencies[k])
                if l == "sini":
                    pattern = self.PatternHelper.getSineIntegerLayout(470e-9, 654e-9, 1e-9, 0.5, 0.5, np.pi, frequencies[k])
                if l == "cosini":
                    pattern = self.PatternHelper.getSineIntegerLayout(470e-9, 654e-9, 1e-9, 0.5, 0.5, -np.pi/2, frequencies[k])
                plotWave = []
                plotInt = []
                #loop over number of iteration steps
                for i in range(23):
                    #loop over number of wavelength steps
                    for j in range(8):
                        #set laser output
                        self.laser.setLaserOutputPM(j, int(pattern.wavelengths[i*8+j]), int(pattern.intensities[i*8+j]/corr_func(pattern.wavelengths[i*8+j]*1e-3)))
                        #append wavelength
                        plotWave.append(pattern.wavelengths[i*8+j])
                        #append intensity
                        plotInt.append(int(pattern.intensities[i*8+j]/corr_func(pattern.wavelengths[i*8+j]*1e-3)))
            
                    #get wavelengths of spectrometer
                    current_wave = self.spec.wavelengths()
                    #get current spectrum
                    if i == 0:
                        current_int = self.spec.intensities()
                        img = self.cam.takeImage()
                        
                    else:
                        current_int += self.spec.intensities()
                        img += self.cam.takeImage()
                        
                plt.plot(current_wave,current_int-np.mean(current_int))
                plt.xlim((470,654))
                plt.show()
        
                #save spectrum
                self.stsToCSV(current_int,current_wave,r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\fourth_image\\"+l+ "_" +str(k)+"nm_STS.spec")
                #save image
                self.imageToCSV(img, r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\fourth_image\\"+l + "_" + str(k)+"nm_CEMOS.data")
                #save pattern
                np.savetxt(r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\fourth_image\pattern_" + l + "_" +str(k)+".csv",(np.asarray(plotWave),np.asarray(plotInt)),delimiter=",")
            
        
        #save frequencies
        np.savetxt(frequencies,r"C:\Users\Labor_NLO-Admin\Desktop\GhostLaser\data\fourth_image\frequencies.csv")
        #get stopping time
        t1 = time.time()
        #get total measurement time
        total = t1-t0
        #print total measurement time
        print("time [s]")
        print(total)
        #plot the desired pattern
        
        fig, ax = plt.subplots(1, gridspec_kw={'wspace': 0}, facecolor='w',figsize=(11,8))
        ax.scatter(plotWave,plotInt, label='Spectrum')
        

    """
    --------------------------------------------------------------------------------------------------------
    --------------------------------Take Curve----------------------------------------------------------
    --------------------------------------------------------------------------------------------------------
    """
    
    @asyncSlot()
    async def startTakeCurve(self):
        startLength = 500
        endLength = 700
        step = 10
        self.spec.integration_time_micros(1000)
        self.cam.setExposure(0)
        self.cam.setTimeOut(100)
        
        
        for i in range(startLength, endLength+step, step):
            self.laser.setLaserOutputPM(0, i*1000, int(10000*self.getCalibratedIntensityFactor(i,True)))
            img = self.cam.takeImage()
            current_wave = self.spec.wavelengths()
            current_int = self.spec.intensities()
            self.imageToCSV(img, "powerCurve/"+str(i)+"nm_CEMOS.data")
            self.stsToCSV(current_int,current_wave,"powerCurve/"+str(i)+"nm_STS.spec")
            schritteGes = (endLength-startLength)/step
            momSchritt = (i - startLength)/step
            percent = momSchritt/schritteGes
            print(percent)
            self.progressCurve.setValue(percent*100)
        
    @asyncSlot()
    async def Take2DCurve(self):
        print("starting 2d curve")
        startLength = 470
        endLength = 660
        step = 5
        self.spec.integration_time_micros(4000)
        self.cam.setExposure(60)
        self.cam.setTimeOut(100)
        startPromille = 0
        endPromille = 1000
        stepPromille = 100
        # Das sind insgesamt 5400 dateien! manuell geht da nichts mehr!
        #zum handeln einfach die zwei for loops kopieren und die dateien einlesen von den Pfaden
        # geschaetzte gesamt groesse fuer die matrix wird 8.2GB !!!
        #daher bitte nicht mehr stuezstellen aufnehmen und wenn dann interpolieren
        for i in range(startLength, endLength+step, step):
            print("wavelength")
            print(i)
            for j in range(startPromille, endPromille+stepPromille, stepPromille):
                print("promille")
                print(j)
            #for i in range(startLength, endLength+step, step):
                self.laser.setLaserOutputPM(0, i*1000, j)
                img = self.cam.takeImage()
                current_wave = self.spec.wavelengths()
                current_int = self.spec.intensities()
                self.imageToCSV(img, "powerCurve_23_09/Pro"+str(j)+"/"+str(i)+"nm_CEMOS.data")
                self.stsToCSV(current_int,current_wave,"powerCurve_23_09/Pro"+str(j)+"/"+str(i)+"nm_STS.spec")
            schritteGes = (endLength-startLength)
            momSchritt = (i - startLength)
            percent = int(momSchritt/schritteGes*100)
            print(percent)
            self.progressCurve.setValue(percent)
        
        
    @asyncSlot()
    async def camerShot(self):
        
        
            toterPixelKoords = (203, 744)
            toterPixel2 = (255, 302)
            now = round(time.time())
            self.cam.setExposure(10)
            img = self.cam.takeImage()
            current_wave = self.spec.wavelengths()
            current_int = self.spec.intensities()
            self.imageToCSV(img, "imgs/"+str(now)+"_CEMOS.data")
            self.stsToCSV(current_int,current_wave,"imgs/"+str(now)+"_SPEC.spec")
            #ja ich lese das von der festplatte gespeicherte ein dann sieht das 1 zu 1 aus wie spaeter auch
            dataImg = pd.read_csv("imgs/"+str(now)+"_CEMOS.data", sep=',')        
            dataImg.drop(dataImg.columns[0], axis=1, inplace=True)
            D = dataImg.to_numpy()
            D = np.transpose(D)
            D[toterPixelKoords] = 0
            D[toterPixel2] = 0

            self.currentPlottedData = D
            self.currentPlottedMaximum = max(map(max, D))
            self.labelCounts.setText(str(self.currentPlottedMaximum))
            self.PlotWidget.canvas.ax.imshow(D, interpolation='none')
            self.PlotWidget.canvas.draw()
            
            self.PlotWidgetSpec.canvas.ax.cla()
            self.PlotWidgetSpec.canvas.draw()
            self.PlotWidgetSpec.canvas.ax.plot(current_wave,current_int)
            self.PlotWidgetSpec.canvas.draw()
            #maxIndex finden fuer die roi
            self.maxIndex = np.unravel_index(D.argmax(), D.shape)
            print(self.maxIndex)
            
            
    @asyncSlot()
    async def fitCurrentShot(self):
        def twoD_Gaussian( xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
            x,y = xdata_tuple
            xo = float(xo)
            yo = float(yo)    
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
            g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                    + c*((y-yo)**2)))
            return g.ravel()
        arr = self.currentPlottedData
    
        x = np.linspace(0, len(arr[0,:]), len(arr[0,:]))
        y = np.linspace(0, len(arr[:,0]), len(arr[:,0]))
        x, y = np.meshgrid(x, y)
    
        maximum = self.currentPlottedMaximum
        initial_guess = (maximum, self.maxIndex[0], self.maxIndex[1], 20, 40, 0, 10)
    
        data_noisy = arr.ravel()
        popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess,maxfev = 10600)#5600
        data_fitted = twoD_Gaussian((x, y), *popt)
        
    
        self.PlotWidget.canvas.ax.cla()
        self.PlotWidget.canvas.draw()    
        self.PlotWidget.canvas.ax.imshow(data_noisy.reshape(len(arr[:,0]), len(arr[0,:])))
        self.PlotWidget.canvas.ax.contour(x, y, data_fitted.reshape(len(arr[:,0]), len(arr[0,:])), 8, colors='w')
        self.PlotWidget.canvas.draw()
        
    @asyncSlot()
    async def roiCurrentShot(self):
        D = self.currentPlottedData
        
        maxCenteredTupel = self.maxIndex
        (xCenter, yCenter) = maxCenteredTupel
        roied = D[xCenter-100:xCenter+100, yCenter-100:yCenter+100]
        self.PlotWidget.canvas.ax.cla()
        self.PlotWidget.canvas.draw()  
        self.currentPlottedData = roied
        self.currentPlottedMaximum = max(map(max, roied))
        self.labelCounts.setText(str(self.currentPlottedMaximum))
        self.PlotWidget.canvas.ax.imshow(roied, interpolation='none')
        self.PlotWidget.canvas.draw()
        
        #neues max enter setzen
        self.maxIndex = np.unravel_index(roied.argmax(), roied.shape)
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    loop = QEventLoop(app)
    asyncio.set_event_loop(loop)
    win = Window()
    win.show()
    
    #with loop:
        #loop.run_forever()

    sys.exit(app.exec())