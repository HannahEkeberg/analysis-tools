from scipy import constants
import periodictable as pt


class NuclearProperties:

    def __init__(self, a):
        self.a = a

    def bindingEnergy(self, Z, N, atom, element_name): # eg. (1, 1, 2H, hydrogen) or (2, 4, 4He, helium)
        atomic_mass_constant = constants.physical_constants["atomic mass constant energy equivalent in MeV"][0]
        # Convert kg to MeV/c2 or atomic mass unit u: [kg] = J/c2 --> 1J = 6242x10^18 MeV/c2 or u
        joule_atomic_mass_constant = constants.physical_constants["joule-atomic mass unit relationship"][0]
        m_p = constants.m_p * joule_atomic_mass_constant # mass proton
        m_n = constants.m_n * joule_atomic_mass_constant # mass neutron
        # element_class = globals()[element_name]
        # element = element_class.isotopes[N]  # Atomic mass 
        # print(elements.name(element_name))
        # atomic_mass = element_class.xxx
        # print(atomic_mass)

        g = pt.H
        isotopes = g.isotopes
        # print(g.mass)
        # pt.mass #.mass('2H')

        g.isotopes('2-H')
        # isotope_197 = g(2).atomic_mass



NuclearProperties("Hi").bindingEnergy(1,2,'2H', 'hydrogen')