import sys
import numpy as np 
import curie as ci 
import pandas as pd
import matplotlib.pyplot as plt
import os

sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
pathToCalibrationSpectra = os.getcwd()

class Calibration:

    def __init__(self, spectra, detectorCalibrationName):
        self.pathToCalibrationSpectrum = ('spectra/calibration' + '/')

        # self.spec_Ba133 = spec_Ba133; self.spec_Cs137 = spec_Cs137; self.spec_Eu152 = spec_Eu152#; self.spec_Am141 = spec_Am141
        self.spectra = spectra
        self.detectorCalibrationName = detectorCalibrationName
        calibrationSpectra = []
        sources = []
        for isotope, spec in self.spectra.items():
            calibrationSpectra.append(self.spectrum(isotope, spec))
            isotope_callable = getattr(self, self.iso_to_method(isotope))
            sources.append(isotope_callable())

        cb = ci.Calibration()
        cb.calibrate(calibrationSpectra, sources=sources)
        cb.plot(show=True)

        cb.saveas(os.getcwd() + '/generatedfiles/calibrationfiles/' + self.detectorCalibrationName + '.json')
        cb.plot(show=True, saveas = './generatedfiles/calibrationfiles/figures/' + self.detectorCalibrationName +'.pdf')

    def spectrum(self, isotope, spec): # isotope = '133BA'
        if type(spec) is list:
            sp = ci.Spectrum(self.pathToCalibrationSpectrum + spec[0]) + ci.Spectrum(self.pathToCalibrationSpectrum + spec[1])
        else:
            sp = ci.Spectrum(self.pathToCalibrationSpectrum + spec)
        sp.isotopes = [isotope]
        return sp
    
    def iso_to_method(self,name: str) -> str:
        digits = ''.join(ch for ch in name if ch.isdigit())
        letters = ''.join(ch for ch in name if ch.isalpha())
        return letters.capitalize() + digits

    def Ba133(self):
        return {'isotope':'133BA', 'A0':4.08480e4, 'ref_date':'01/11/2011 12:00:00'} # todo fix, from uio

    def Co60(self):
        return {'isotope':'60CO', 'A0':3.3859E4, 'ref_date':'03/01/2019 12:00:00'} # todo fix, from uio

#det 80:
spectra_det80_0cm = {'60CO': 'AB18052026_60Co_det80_0cm.Spe', '133BA': 'AA18052026_133Ba_det80_0cm.Spe'}
Calibration(spectra_det80_0cm, 'calibration_det80_0cm')

# det 40
# spectra_det50_3_6cm = {'133BA': 'AA18052025_Ba133_det50_3_6cm.Spe'}
# Calibration('spectra_det50_3_6cm', spectra_det50_3_6cm, 'spectra_det50_3_6cm')

# spectra_det2_10 = {'137CS': 'BA20251017_Det2_Cs137_10cm.Spe', '133BA': 'BV20251017_Det2_Ba133_10cm.Spe', '152EU': 'CL20251020_Det2_Eu152_10cm.Spe'}
# Calibration('Det2_10', spectra_det2_10, 'calibration_det2_10')