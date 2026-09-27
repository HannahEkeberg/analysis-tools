import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import os
import sys

# sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
# from nuclearanalysistools.Tendl import *
import curie as ci

class AnalyzeSpectrum_dy:
    def __init__(self, detector, calibrationFile, peakDataFolderName=None):
        self.calibrationFile = calibrationFile
        self.detector = detector
        self.peak_data_folder_name = '/' + peakDataFolderName + '/' if peakDataFolderName != None else '/data/'
        self.pathToSpectra = os.getcwd() + '/spectra/'
        self.pathToCalibrationFiles = os.getcwd() + '/generatedfiles/calibrationfiles/'
        self.pathToPeakData = os.getcwd() + '/generatedfiles/peakdata/'
        
    def analyze(self, spectrumFileName, listOfIsotopes, x_rays=False):
        # spectrumFileName = spectrumFileName
        if spectrumFileName.endswith(".Spe"):
            spectrumName = spectrumFileName.replace(".Spe", "")
        else:
            spectrumName = spectrumFileName
        try:
            cb = ci.Calibration(self.pathToCalibrationFiles + self.calibrationFile)
            sp = ci.Spectrum(self.pathToSpectra + self.detector + '/' + spectrumFileName)
            # sp.cb = cb
            sp.isotopes = listOfIsotopes
            if x_rays:
                sp.fit_config = {'xrays':True, 'E_min':20}
            else:
                sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
            sp.plot()
            sp.saveas(self.pathToPeakData + self.peak_data_folder_name + spectrumName + '_peak_data.csv', replace=True)
            sp.saveas(self.pathToPeakData + '/figures/' + spectrumName + '_peak_data.pdf', replace=True) 
        except AttributeError:
            print("Most likely no peaks...")

    def analyze_jobs(self, listOfJobSpectra, listOfIsotopes, peakSummaryFilename=None):
        if peakSummaryFilename == None:
            peakSummaryFilename = listOfJobSpectra[0].replace("000.Spe", "job")
        cb = ci.Calibration(self.pathToCalibrationFiles + self.calibrationFile)
        spectra = []
        for job in listOfJobSpectra:
            sp = ci.Spectrum(self.pathToSpectra + self.detector + '/' + job)
            spectra.append(sp)

        summed_spectrum = spectra[0]
        for i, spec in enumerate(spectra[1:], start=1):
            summed_spectrum += spec
            print(f"Added spectrum {i+1}/{len(spectra)}")

        # summed_spectrum.cb = cb
        summed_spectrum.isotopes = listOfIsotopes
        # summed_spectrum.fit_config = {'SNR_min':3.5}
        summed_spectrum.fit_config = {'xrays':False, 'E_min':20}
        summed_spectrum.plot()
        try:
            pass
            summed_spectrum.saveas('generatedFiles/peakData/figures/' + peakSummaryFilename + '_peak_data.pdf')
            summed_spectrum.saveas('generatedFiles/peakData/data/' + peakSummaryFilename + '_peak_data.csv', replace=False)
        except AttributeError:
            print("Most likely no peaks...")
