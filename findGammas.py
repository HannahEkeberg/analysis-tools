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
        if not data.empty:
            print(data)
        else:
            print(r'No matching decay gammas for: ' + str(gammaLine) + ' (+/- ' + str(gammaLineTolerance) + ') keV')
        return data

    def findAllGammas(self, minIntensity=None, xrays=False):
        gammas = [] # ('63Zn', thalf, energy, intensity)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            gammaLines = isotope.gammas(I_lim=minIntensity, xrays=xrays)#, dE_511=1.0)
            if not gammaLines.empty:
                for i in range(len(gammaLines)):
                    energy = gammaLines.at[i, 'energy']; intensity=gammaLines.at[i, 'intensity']
                    gammas.append([iso, energy, intensity, half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(
            gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Energy')
        print(data)
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
        print(data)
        return data

    def orderIsotopesByHalfLife(self):
        isotopeInfo = [] # isotope, t half (s,m,h,d)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            isotopeInfo.append([iso, half_life, half_life/60, half_life/(60*60), half_life/(60*60*24)])
        data = pd.DataFrame(
            isotopeInfo, columns=['Isotope', 'Half life (s)', 'Half life (m)', 'half life (h)', 'half life (d)']).sort_values('Half life (s)')
        print(data)
        return data

    def saveCsv(self, dataframe, filename, directory=None):
        if directory==None:
            filename = os.getcwd() + '/' + filename + '.csv'
        else:
            filename = os.getcwd() + '/' + directory + '/' + filename + '.csv'
        dataframe.to_csv(filename, index=False)





#EXAMPLE

listOfPossibleProductsNatCu= ['61ZN', '62ZN', '63ZN', '65ZN', '60CU', '61CU', 
                            '62CU', '64CU', '56NI', '57NI', '63NI', '65NI',
                            '55CO', '56CO', '57CO', 'CO-58m', 'CO-58g', '60CO',
                            '61CO', '55FE', '59FE']
## Can also add typical isotopes causing background radiation

ag = AnalyzeGammas(listOfPossibleProductsNatCu)
# ag.matchByGamma(2614.0, gammaLineTolerance=2.5, minIntensity=None)
allGammas = ag.findAllGammas(minIntensity=0.1)
ag.saveCsv(allGammas, 'gammas_natCu', 'csv')
# ag.findGammasSpecificIsotope('63ZN')
# ag.orderIsotopesByHalfLife()


