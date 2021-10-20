# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 14:44:05 2021

@author: Jonas Brandhoff
"""
from NKTP_DLL import *

class Laser:
    def __init__(self):
        self.port = 'COM5'
        self.open  = False
        self.mainLaserAd = 15
        self.RFAd = 16
        self.doPrint = False
    def openConnection(self):
        openResult = openPorts(self.port, 0, 0)
        if self.doPrint: print('Opening the comport:', PortResultTypes(openResult))
        self.open = True
        
    def closeConnection(self):
        closeResult = closePorts('COM5')
        self.open = False
    def EmissionOn(self):
        if not self.open:
            if self.doPrint: print("please connect first")
            return
        result = registerWriteU8(self.port, self.mainLaserAd, 0x30, 2, -1)
        if self.doPrint: print('Setting emission ON - Extreme:', RegisterResultTypes(result))
        
    def EmissionOff(self):
        registerWriteU8(self.port, self.mainLaserAd, 0x30, 0x00, -1)
        
    
    
    def resetInterLock(self):
        result = registerWriteU8(self.port, self.mainLaserAd,0x32,1,-1)
        if self.doPrint: print('Resetting interlock', RegisterResultTypes(result))
    
    
    
    
    def RFOn(self):
        result = registerWriteU8(self.port, self.RFAd, 0x30, 1, -1)
        if self.doPrint: print('changeing amplitude', RegisterResultTypes(result))
    def RFOff(self):
        result = registerWriteU8(self.port, self.RFAd, 0x30, 0, -1)
        if self.doPrint: print('changeing amplitude', RegisterResultTypes(result))
    
    
    
    def activateFSK(self):
        result = registerWriteU8(self.port, self.RFAd, 0x3B, 1, -1)
        if self.doPrint: print('activating FSK', RegisterResultTypes(result))
    
    def deactivateFSK(self):
        result = registerWriteU8(self.port, self.RFAd, 0x3B, 0, -1)
        if self.doPrint: print('deactivating FSK', RegisterResultTypes(result))
    
    
    def changeAmplitude(self, index, valueInPercent):
        Adress = 0xB0
        Adress = Adress + index
        result = registerWriteU16(self.port, self.RFAd, Adress, valueInPercent*10, 0)
        if self.doPrint: print('changeing amplitude', RegisterResultTypes(result))

    def changeAmplitudePER(self, index, valueInPermill):
        Adress = 0xB0
        Adress = Adress + index
        result = registerWriteU16(self.port, self.RFAd, Adress, valueInPermill, 0)
        if self.doPrint: print('changeing amplitude', RegisterResultTypes(result))



    def changeWavelength(self, index, valueInNm):
        Adress = 0x90
        Adress = Adress + index
        if self.doPrint: print('Adress is:'+str(Adress))
        result = registerWriteU32(self.port, self.RFAd, Adress, valueInNm*1000, 0)
        if self.doPrint: print('changeing wavelength', RegisterResultTypes(result))
        
    def changeWavelengthPM(self, index, valueInPm):
        Adress = 0x90
        Adress = Adress + index
        if self.doPrint: print('Adress is:'+str(Adress))
        result = registerWriteU32(self.port, self.RFAd, Adress, valueInPm, 0)
        if self.doPrint: print('changeing wavelength', RegisterResultTypes(result))
        

        
    def setLaserOutput(self, index, waveLength, amplitude):
        self.changeWavelength(index, waveLength)
        self.changeAmplitude(index, amplitude)
    
    def setLaserOutputPM(self, index, waveLengthPM, amplitudePER):
        self.changeWavelengthPM(index, waveLengthPM)
        self.changeAmplitudePER(index, amplitudePER)