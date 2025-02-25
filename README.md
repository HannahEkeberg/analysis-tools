# analysis-tools

Tools for plotting reaction modeling codes


findGammas:
	Tool to:
- match by specific gamma.
- find all gammas in pre-defined list of isotopes
- find all gammas for specified isotope
- order all pre-defined isotopes by half life
- save csv. All methods above returns data frames which can be saved as csv.
	
countingTimes:
	Tool to: 
- get counting statistics. Eg. How long to count an isotopes with a specific energy and intensity (based on efficiency, A0, how many counts you want). Takes in a dictionary of isotopes with energies of interest and corresponding intensity. Returns a dataframe. 
- save csv



