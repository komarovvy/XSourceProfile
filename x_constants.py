#TODO classes of absorption and radiations
# characteristic wavelengthes, [Ang]
WAVELENGTH = {
    'MoKaw': 0.71075,
    'MoKa1': 0.70932, 
    'MoKa2': 0.71361,
    'CuKaw': 1.54187,
    'CuKa1': 1.54059, 
    'CuKa2': 1.54443,
    }

# name and density, [g/cm**3]
MATERIALS = {
    'Fe': ('alfa-Fe, room T', 7.86)
    }

# linear absorption coefficients for the spec. material and w.length, [mm**-1]
# data from NIST, https:\\.... #TODO find the URL!
#TODO find more convenient way to represent absorption coeffs for materials/wl
PHASE_TO_MU = {
    'Fe, MoKaw': (36.173*2 + 36.781*1)/3 * 7.86 * 0.1,
    'Fe, MoKa1': 36.173*7.86*0.1,
    'Fe, MoKa2': 36.781*7.86*0.1,
    'Fe, CuKaw': (300.06*2 + 302.12*1)/3 * 7.86 * 0.1,
    'Fe, CuKa1': 300.06*7.86*0.1,
    'Fe, CuKa2': 302.12*7.86*0.1,
    }