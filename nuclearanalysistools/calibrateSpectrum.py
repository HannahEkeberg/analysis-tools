import curie as ci
import pandas as pd

class Calibrate:

    def __init__(self, pathToSpectra, source, calibrationSpectrumNames):
        self.cb = ci.Calibration()
        self.pathToSpectra = pathToSpectra
        self.sourcs = sources #  [ {'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'}, {'isotope':'152EU', 'A0':370000, 'ref_date': '11/01/1984 12:00:00'}]
        self.calibrationSpectrumNames = calibrationSpectrumNames # []
        

    def plot(self):
        pass

    def saveEfficiency(self):
        pass