# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 16:54:47 2021

@author: Paul Herrmann
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.interpolate as itp
from numpy import fft

#first step: preprocessing: process data from camera, either roied or not

def get_file_names(input_path):
    #get list of all files in input folder
    files = os.listdir(input_path)
    spectra_files = []
    image_files = []
    wavelengths = []
    pattern_files = []
    for i in range(len(files)):
        #append all image files to image files list
        
        #this means that it is an image file from the camera
        if files[i].endswith("a") == True:
            image_files.append(input_path + "/" +files[i])
            #cosine input pattern
            if files[i].startswith("c") == True:
                wavelengths.append(float(files[i][7]))
            #sine input pattern
            if files[i].startswith("s") == True:
                wavelengths.append(float(files[i][5]))
        #this is a csv file containing the input pattern
        if files[i].startswith("pattern") == True:   
            pattern_files.append(input_path + "/" + files[i])      
        #append all spectrum files to spectrum files list
        if files[i].endswith("c") == True:
            spectra_files.append(input_path + "/" + files[i])
        else:
            continue
    return spectra_files, image_files,pattern_files, wavelengths

def preprocess_images(spectra_files,image_files,output_path,wavelengths,roied=False,plotting=False):
    print("Step 1:Preprocessing")
    bucket = np.zeros(len(spectra_files))
    #get the wavelengths as a list of floats
    for i in range(len(spectra_files)):
        dataImg = pd.read_csv(image_files[i], sep=',')
        dataImg.drop(dataImg.columns[0], axis=1, inplace=True)
        D = dataImg.to_numpy()
        D = np.transpose(D)
        if roied == True:
            bucket[i] = np.sum(D[300:600,300:700])
        else :
            bucket[i] = np.sum(D)
    #save bucket intensity to file output_path
    np.savetxt(output_path,(np.array(wavelengths),bucket),delimiter=",",header="wavelength [nm], bucket intensity")
    if plotting == True:
        plt.plot(wavelengths, bucket)
        plt.show()
    print("Preprocessing finished\n")

def ghost_imaging(spectra_files, bucket_file,pattern_files,save_path,method="TGI",plotting=False,pattern=False):
    print("\t Starting ghost spectrum recovery with " + method + " method")
    methods = ["TGI", "DGI","NGI","PIGI"]
    if method not in methods:
        print("error, please select a correct ghost imaging method!")
    #load bucket
    wavelengths, bucket = np.loadtxt(bucket_file,delimiter=",")
    #get wavelengths for spectra
    if pattern == False:
        data = pd.read_csv(spectra_files[0], sep=',')
        data.drop(data.columns[0], axis=1, inplace=True)
        spectra_wavelengths  = data.to_numpy()[308:708,1]
    else:
        spectra_wavelengths, spectrum = np.loadtxt(pattern_files[0],delimiter=",")
    #initialize arrays for avg input and reference spectra
    avg_input_spectrum = np.zeros(len(spectra_wavelengths))
    avg_reference = 0
    recovered_spectrum = np.zeros(len(spectra_wavelengths))
    integrated_spectrum = np.zeros(len(bucket))
    avg_bucket = 0
    
    #get the avg values
    for i in range(len(bucket)):
        if pattern == False:
            data = pd.read_csv(spectra_files[i], sep=',')
            data.drop(data.columns[0], axis=1, inplace=True)
            spectrum = data.to_numpy()[308:708,0]
        else:
            ws,spectrum = np.loadtxt(pattern_files[i],delimiter=",")
            spectrum /= 1000
            if i == 0:
                plt.plot(ws,spectrum)
                plt.show()
        avg_input_spectrum += spectrum 
        avg_bucket += bucket[i]
        avg_reference += np.trapz(spectrum)
        integrated_spectrum[i] = np.trapz(spectrum)
    avg_input_spectrum /= len(bucket)
    avg_bucket /= len(bucket)
    avg_reference /= len(bucket)
    
    #ghost imaging algorithm
    if method == "PIGI":
        A = np.zeros((len(bucket),len(spectra_wavelengths)))
    for i in range(len(bucket)):
        #load spectrum
        if pattern == False:
            data = pd.read_csv(spectra_files[i], sep=',')
            data.drop(data.columns[0], axis=1, inplace=True)
            spectrum = data.to_numpy()[308:708,0]
        else:
            ws,spectrum = np.loadtxt(pattern_files[i],delimiter=",")
        
        if method == "PIGI":
            A[i] = spectrum - avg_input_spectrum
        
        if method == "TGI":
            recovered_spectrum += (spectrum - avg_input_spectrum)*(bucket[i]-avg_bucket)
        if method == "DGI":
            recovered_spectrum += (spectrum - avg_input_spectrum)*(bucket[i]-avg_bucket/avg_reference*np.trapz(spectrum))
        if method == "NGI":
            recovered_spectrum += (spectrum - avg_input_spectrum)*(bucket[i]/np.trapz(spectrum)-avg_bucket/avg_reference)
    recovered_spectrum /= len(bucket)
    if method == "PIGI":
        
        recovered_spectrum = np.matmul(np.linalg.pinv(A),bucket-avg_bucket/avg_reference*integrated_spectrum)
    
    
    
    if plotting == True:
        plt.plot(spectra_wavelengths,recovered_spectrum)
        plt.xlabel("wavelength [nm]")
        plt.ylabel("intensity [a.u.]")
        plt.title(method)
        plt.show()
        
    save_file = save_path + "/" + method + ".gspec"
    np.savetxt(save_file,(spectra_wavelengths,recovered_spectrum),delimiter=",",header="wavelengths [nm], intensity")
    print("\t Ghost spectrum recovery finished")
    return save_file


def OCT(spectrum_file,method,zero_padding=False,plotting=False):
    print("\t Starting OCT of " + method +" method\n")
    #load data
    wavelengths, spectrum = np.loadtxt(spectrum_file,delimiter=",")
    spectrum *= np.hamming(len(spectrum))
    #define k
    k = np.pi/wavelengths
    #interpolate
    fct = itp.interp1d(k,spectrum)
    #define equidistant k    
    k_eq = np.linspace(k[-1],k[0],len(k))
    #interpolate on equidistant grid
    spectrum_inter = fct(k_eq)
    if zero_padding == True:
        spectrum_inter = np.zeros(3*len(k_eq))
        spectrum_inter[:len(k_eq)] = fct(k_eq)
        k_eq = np.linspace(k[-1],k[0]+2*(k[0]-k[-1]),3*len(k))
        
    x = np.fft.fftfreq(len(k_eq),k_eq[1]-k_eq[0])*1e-9
    trans = np.abs(np.fft.fft(spectrum_inter-np.mean(spectrum_inter)))
    
    
    np.savetxt(spectrum_file + "x", (x,trans),delimiter=",",header="x [m], intensity")
    if plotting == True:
        plt.title(method)
        plt.plot(x,trans)
        plt.xlabel("x [m]")
        plt.ylabel("OCT signal [a.u.]")
        #plt.xlim((0,25e-6))
        plt.show()
        
def direct_pattern(spectra_files, bucket_file, output_folder):
    wavelengths, bucket = np.loadtxt(bucket_file,delimiter=",")
    l = len(bucket)/4
    x = np.linspace(0,l-1,l)
    #get the cosine coefficients
    a = bucket[0:l] - bucket[l:2*l]
    #get the sine coefficients
    b = bucket[2*l:3*l] - bucket[3*l:4*l]
    #absolute value
    transformed = np.sqrt(a**2+b**2)
    #plot the result
    plt.plot(x,transformed)
    plt.show()
    #save the result
    np.savetxt(output_folder + "/direct.csv",(x,transformed),delimiter=",")
    
    
def ghost_complete(input_folder,output_folder,preprocessing=False,pattern=False,direct=False):
    #input_folder, string, is the path of the data that has to be processed
    #output_folder, string, is the path where the results should be saved
    #preprocessing, boolean, defines if the camera pics have to be integrated to get the bucket intensity, can be waived if bucket file exists
    #pattern, boolean, if True the pattern from the pattern file instead of the measured input spectrum is taken for ghost spectrum recovery, i.e. quasi-computational --> computational
    
    #available methods for ghost spectrum recovery
    methods = ["TGI","DGI","NGI","PIGI"]
    #define a file for the bucket intensity
    bucket_file = input_folder + "/bucket.csv"
    #get all necessary file lists
    spectra_files, image_files,pattern_files,wavelengths = get_file_names(input_folder)
    #preprocess the data
    if preprocessing == True:
        preprocess_images(spectra_files,image_files,bucket_file,wavelengths)
    #ghost imaging
    print("Ghost OCT")
    #loop over all possible methods
    for i in range(len(methods)):
        #perform ghost spectrum recovery
        save_file = ghost_imaging(spectra_files,bucket_file,pattern_files,output_folder,method=methods[i],plotting=True,pattern=pattern)
        #perform FFT to get OCT result
        OCT(save_file,methods[i],zero_padding=True,plotting=True)
    if direct == True:
        direct_pattern(spectra_files, bucket_file,output_folder)
    print("Finished data evaluation")
    
#ghost_complete(r"C:\Users\Paul\Desktop\third_image",r"C:\Users\Paul\Desktop\third_result",preprocessing=False,pattern=False,direct=False)

inde, intensity,wavelengths, = np.loadtxt(r"C:\Users\Paul\Desktop\third_image\cosini_1.18877001551221e-14nm_STS.spec",skiprows=1,delimiter=",",unpack=True)
plt.plot(wavelengths,intensity)
plt.show()