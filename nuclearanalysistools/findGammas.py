import curie as ci
import pandas as pd
import os

class AnalyzeGammas:

    def __init__(self, isotopes=[]):
        self.isotopes = isotopes # a list of possible nuclei that can be matched in spectra

    def matchByGamma(self, gammaLine, gammaLineTolerance=0.5, minIntensity=None, xrays=False):
        Elim = (gammaLine - gammaLineTolerance, gammaLine + gammaLineTolerance )
        possibleGammas = [] # ('63Zn', thalf, energy, intensity)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            matchedGamma_dataframe = isotope.gammas(I_lim=minIntensity, xrays=xrays, E_lim=Elim)#, dE_511=1.0)
            half_life = isotope.half_life()
            if not matchedGamma_dataframe.empty:
                for i in range(len(matchedGamma_dataframe)):
                    energy = matchedGamma_dataframe.at[i, 'energy']; intensity=matchedGamma_dataframe.at[i, 'intensity']
                    possibleGammas.append([iso, str(energy), str(intensity), half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(possibleGammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Half life (s)')
        if data.empty:
            raise Exception('No matching decay gammas for: ' + str(gammaLine) + ' (+/- ' + str(gammaLineTolerance) + ') keV')
        else:
            return data

    def findAllGammas(self, minIntensity=None, xrays=False):
        gammas = [] # ('63Zn', thalf, energy, intensity)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            gammaLines = isotope.gammas(I_lim=minIntensity, xrays=xrays)#, dE_511=1.0)
            if not gammaLines.empty:
                for i in range(len(gammaLines)):
                    energy = gammaLines.at[i, 'energy']; intensity=gammaLines.at[i, 'intensity']; unc_intensity=gammaLines.at[i, 'unc_intensity']
                    gammas.append([iso, energy, intensity,unc_intensity,half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(
            gammas, columns=['Isotope', 'Energy', 'Intensity', 'unc_intensity', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Energy')
        return data

    def findGammasSpecificIsotope(self, iso, minIntensity=None, xrays=None):
        isotope = ci.Isotope(iso)
        half_life = isotope.half_life()
        gammas = []
        gammaLines = isotope.gammas(I_lim=minIntensity, xrays=xrays)#, dE_511=1.0)
        if not gammaLines.empty:
            for i in range(len(gammaLines)):
                energy = gammaLines.at[i, 'energy']; intensity=gammaLines.at[i, 'intensity']
                gammas.append([iso, energy, intensity, half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(
            gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Intensity', ascending=False)
        return data

    def orderIsotopesByHalfLife(self):
        isotopeInfo = [] # isotope, t half (s,m,h,d)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            isotopeInfo.append([iso, half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(
            isotopeInfo, columns=['Isotope', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Half life (s)')
        return data

    def saveCsv(self, dataframe, filename, directory=None):
        if directory==None:
            filename = os.getcwd() + '/' + filename + '.csv'
        else:
            filename = os.getcwd() + '/' + directory + '/' + filename + '.csv'
        dataframe.to_csv(filename, index=False)


# AnalyzeGammas([]).findGammasSpecificIsotope('IR-192g')


#EXAMPLE

# listOfPossibleProductsNatCu= ['61ZN', '62ZN', '63ZN', '65ZN', '60CU', '61CU', 
#                             '62CU', '64CU', '56NI', '57NI', '63NI', '65NI',
#                             '55CO', '56CO', '57CO', 'CO-58m', 'CO-58g', '60CO',
#                             '61CO', '55FE', '59FE']

# listOfPossibleProductsNatIr= ['PT-193m', '191PT', '189PT', '188PT', '187PT', '186PT', 
#                             'IR-194m', 'IR-194g', 'IR-192m2', 'IR-192m1', 'IR-192g', 'IR-190m2',
#                             'IR-190m1', 'IR-190g', '189IR', '188IR', 'IR-187g', '193OS',
#                             '191OS']

# listOfPossibleProductsNatFe= ['55CO', '56CO', '57CO', 'CO-58m', 'CO-58g', '60CO',
#                             '61CO', '55FE', '59FE', 'FE-53m', 'FE-53g', '52FE',
#                             '56MN', '54MN', '53MN', 'MN-52m', 'MN-52g', '51MN',
#                             '56CR', '55CR', '51CR', '49CR', '48CR', '52V', '49V', '48V', '47V']


# listOfPossibleProductsNatTi = ['52V', '49V', '48V', '47V','51TI', '45TI','44TI','49SC',
#                                 '48SC', '47SC', '46SC', '44SC', '43SC',
#                                 '49CA', '47CA', '45CA', '45K', '44K', '43K', '42K', '40K', '24NA']



# from Tendl import *
# tendl = Tendl({"Fe54": 0.05845, "Fe56": 0.91754, "Fe57": 0.02119, "Fe58": 0.00282})
# # tendl.plotTendl23Unique(productZ='26', productA='53', isomerLevel = None, color=None, lineStyle=None, label='53Fe')
# # tendl.plotTendl23Unique(productZ='25', productA='52', isomerLevel = '01', color='red', lineStyle=None, label='52mMn')
# plt.legend()
# plt.show()


# listOfCalibrationSources = ['137CS', '133BA', '152EU']
## Can also add typical isotopes causing background radiation
# ag = AnalyzeGammas(listOfPossibleProductsNatCu)
# ag = AnalyzeGammas(listOfCalibrationSources)

# ag = AnalyzeGammas(listOfPossibleProductsNatFe)

# ag.orderIsotopesByHalfLife()
# ag.findGammasSpecificIsotope('56CO', minIntensity=0.001, xrays=True)


# ag.matchByGamma(2614.0, gammaLineTolerance=2.5, minIntensity=None)
# allGammas = ag.findAllGammas(minIntensity=0.1, xrays=True)
# ag.matchByGamma(377.900, gammaLineTolerance=2.5, minIntensity=None)
# ag.matchByGamma(1434.060, gammaLineTolerance=2.5, minIntensity=None)
# ag.matchByGamma(1727, gammaLineTolerance=2.5, minIntensity=None)
# ag.saveCsv(allGammas, 'gammas_natTi', 'csv')
# ag.findGammasSpecificIsotope('53FE', minIntensity=0.001, xrays=True)
# ag.findGammasSpecificIsotope('MN-52m', minIntensity=0.001, xrays=False)
# ag.orderIsotopesByHalfLife()


