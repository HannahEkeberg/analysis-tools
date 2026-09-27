import sys
# sys.path.append('/opt/homebrew/lib/python3.13/site-packages/curie')
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
            formatted_half_life = self.format_time(half_life)
            if not matchedGamma_dataframe.empty:
                for i in range(len(matchedGamma_dataframe)):
                    energy = matchedGamma_dataframe.at[i, 'energy']; 
                    intensity=matchedGamma_dataframe.at[i, 'intensity']; unc_intensity = matchedGamma_dataframe.at[i, 'unc_intensity']
                    possibleGammas.append([iso, str(energy), str(intensity), str(unc_intensity), formatted_half_life, half_life])
        data = pd.DataFrame(possibleGammas, columns=['Isotope', 'Energy', 'Intensity', 'Unc intensity', 'Half life', 'Half life (s)'])#.sort_values('Intensity')
        if not data.empty:
            print(data)
            return data
        else:
            print('No matching decay gammas for: ' + str(gammaLine) + ' (+/- ' + str(gammaLineTolerance) + ') keV')
            return data

    def findAllGammas(self, minIntensity=None, xrays=False):
        gammas = [] # ('63Zn', thalf, energy, intensity)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            formatted_half_life = self.format_time(half_life)
            gammaLines = isotope.gammas(I_lim=minIntensity, xrays=xrays)#, dE_511=1.0)
            if not gammaLines.empty:
                for i in range(len(gammaLines)):
                    energy = gammaLines.at[i, 'energy']; intensity=gammaLines.at[i, 'intensity']; unc_intensity=gammaLines.at[i, 'unc_intensity']
                    gammas.append([iso, energy, intensity, formatted_half_life, half_life])
        data = pd.DataFrame(
            # gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life', 'Half life (s)']).sort_values(['Energy'], ascending=[True])
            gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life', 'Half life (s)']).sort_values(['Intensity'], ascending=[False])
            # gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life', 'Half life (s)']).sort_values(['Half life (s)', 'Intensity'], ascending=[True, False])
        return data
    
    def matchByAllGammas_decaydb(self, gammaLine, gammaLineTolerance=2.5, minIntensity=None, xrays=False):
        df = ci.Isotope('56NI').get_all_gammas()
        Elim = (gammaLine - gammaLineTolerance, gammaLine + gammaLineTolerance )
        matching_rows = df[df['energy'].between(
            gammaLine - gammaLineTolerance, 
            gammaLine + gammaLineTolerance)]
        matching_rows['half_life'] = matching_rows['isotope'].apply(lambda x: self.format_time(ci.Isotope(x).half_life()))
        print(matching_rows)
        print(matching_rows['isotope'].values)
        return matching_rows


    def findGammasSpecificIsotope(self, iso, minIntensity=0.01, xrays=None):
        isotope = ci.Isotope(iso)

        half_life = isotope.half_life()
        formatted_half_life = self.format_time(half_life)
        gammas = []
        gammaLines = isotope.gammas(I_lim=minIntensity, xrays=xrays, dE_511=1.0)
        if not gammaLines.empty:
            for i in range(len(gammaLines)):
                energy = gammaLines.at[i, 'energy']; intensity=gammaLines.at[i, 'intensity']
                gammas.append([iso, energy, intensity, formatted_half_life])
        data = pd.DataFrame(
            gammas, columns=['Isotope', 'Energy', 'Intensity', 'Half life']).sort_values('Intensity', ascending=False)
        return data

    def orderIsotopesByHalfLife(self):
        isotopeInfo = [] # isotope, t half (s,m,h,d)
        for iso in self.isotopes:
            isotope = ci.Isotope(iso)
            half_life = isotope.half_life()
            formatted_half_life = self.format_time(half_life)
            isotopeInfo.append([iso, formatted_half_life, half_life])
        data = pd.DataFrame(
            isotopeInfo, columns=['Isotope', 'Half life', 'Half life (s)']).sort_values('Half life (s)')
        return data

    def saveCsv(self, dataframe, filename, directory=None):
        if directory==None:
            filename = os.getcwd() + '/' + filename + '.csv'
        else:
            filename = os.getcwd() + '/' + directory + '/' + filename + '.csv'
        dataframe.to_csv(filename, index=False)
    
    def format_time(self, t_seconds):
        if t_seconds < 60:
            return f"{t_seconds:.2f} s"
        elif t_seconds < 3600:
            return f"{t_seconds/60:.2f} min"
        elif t_seconds < 86400:
            return f"{t_seconds/3600:.2f} h"
        elif t_seconds < 31557600:  # ~1 year (365.25 days)
            return f"{t_seconds/86400:.2f} d"
        else:
            return f"{t_seconds/31557600:.2f} y"


