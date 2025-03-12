# analysis-tools

Tools for plotting reaction modeling codes


Check if and where a python package is installed
pip3 show nuclearanalysistools

------------------------------------------------------------------------------------

What makes a python package:
------------------------------------------------------------------------------------
mypackage/
│── mypackage/       # Main package folder (same name as the package)
│   │── __init__.py  # Marks this directory as a package
│   │── module1.py   # A Python module (script)
│   │── module2.py   # Another module
│
│── setup.py         # (Optional) If you want to install it
│── README.md        # (Optional) Documentation
│── requirements.txt # (Optional) Dependencies

------------------------------------------------------------------------------------

Put this into setup.py
from setuptools import setup, find_packages

------------------------------------------------------------------------------------

setup(
    name="mypackage",
    version="0.1",
    packages=find_packages(),
)
------------------------------------------------------------------------------------
In terminal, cd into the package (top level)
pip install .
------------------------------------------------------------------------------------

Pandas tools.

Tendl:

findGammas: Need to provide a list of isotopes that are of interest in curie notation
	Tool to:
- matchByGamma(gammaLine, gammaLineTolerance=0.5, minIntensity=None, xrays=False) --> match by specific gamma, find all mathcing for different nuclei. Must provide list in initialization
- findAllGammas(minIntensity=None, xrays=False) --> find all gammas in pre-defined list of isotopes
- findGammasSpecificIsotope(iso, minIntensity=None, xrays=None) --> find all gammas for specified isotope. Not necessary to pre-define isotopes. Isotope (iso) must have curie notation
- orderIsotopesByHalfLife() --> order all pre-defined isotopes by half life
- saveCsv(dataframe, filename, directory=None) --> save csv. All methods above returns data frames which can be saved as csv.
	
countingTimes: Need to provide foil ('SC'), target {'Sc45': 1.0}, numberOftargetNuclei (ask chat gpt, for instance specify area and thickness of target), beam energy (30 #MeV) and beamParticle ('proton'). 
Assuming: efficiency data from hpge detector in room 131 at 5 cm. assuming beam current 100 nA, irradiation time 2 hrs, desired number of counts is 10 000. 
	Tool to: 
- A0(product, productZ, productA, isomerLevel=None): Estimate A0 as a function of product (eg. SC45), Z, A, isomerlevel. Calculates based on tendl cross sections.
- estimatedA(A0, decay_const, t_delay): Estimate A based on A0, decay constant and time delay since EOB
- N_T( areal_density, molar_mass): Calculate number of target nuclei based on areal density and molar mass of target material. 
- getCountingStatistics(isotopeDict). Takes in isotope dict with {isotope: [[[E1, I1], [E2, I2], ...], Z, A, isomerLevel]}. Gives a dataframe with isotope, half life, time delay, estimated activity, foil, gamma energy, gamma intensity, counting time (min).



