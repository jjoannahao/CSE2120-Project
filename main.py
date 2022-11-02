"""
title: Precipitate Calculator
author: Joanna Hao
date-created: 2022-10-20
"""
# ---------- VARIABLES ---------- #
"""
SOLUBILITY TABLE
column 1 = pos1 ions, neg 1 ions
column 2, 3 = neg 1 ions
column 4, 5 = neg 2 ions
column 6, 7 = neg 1 ions
"""
POS1_IONS = ("H", "Li", "Na", "K", "Rb", "Cs", "Fr",
             "NH4",
             "Ag", "Tl")
POS2_IONS = ("Be", "Mg", "Ca", "Sr", "Ba", "Ra",
             "Mn", "Co", "Ni", "Pd", "Cu", "Zn", "Cd", "Pb", "Po",
             "Pd", "Cd", "Hg2", "Hg", "Po")
POS3_IONS = ("Sc", "Y", "La", "Ac",
             "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
             "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Lr",
             "Cr", "Au", "Bi",
             "Fe", "Ru", "Rh", "Ga", "In", "Sb", "Al")
POS4_IONS = ("Ti", "Zr", "Hf", "Rf",
             "Th", "Pu", "Os", "Ir",
             "Ge", "Sn", "Pt", )
POS5_IONS = ("Pa", "Np",
             "V", "Nb", "Ta", )
POS6_IONS = ("U", "Mo", "W")
POS7_IONS = ("Tc", "Re")
NEG1_IONS = ("F", "Cl", "Br", "I", "At",
             "CH3COO", "C6H5COO", "HCO3", "ClO4", "ClO3", "ClO2", "OCl", "ClO", "CN", "OH", "IO3", "NO3", "NO2", "HOOCCOO", "MnO4", "H2PO4", "HSO4", "HSO3", "HS", "SCN")
NEG2_IONS = ("O", "S", "Se", "Te",
             "C2", "CO3", "CrO4", "Cr2O7", "OOCCOO", "O2", "S2", "HPO4", "SiO3", "SO4", "SO3", "S2O3")
NEG3_IONS = ("N", "P", "As")
NEG4_IONS = ("C", "Si")  # taught that B, C, Si don't tend to form ions

COL2_PRECIPITATES = ("Li", "Mg", "Ca", "Sr", "Ba", "Fe", "Hg2", "Pb")
COL3_PRECIPITATES = ("Cu", "Ag", "Hg2", "Pb", "Tl")  # Cu = +1 charge (less common)
COL4_PRECIPITATES = ("Ca", "Sr", "Ba", "Ag", "Hg2", "Pb", "Ra")
COL5_TO_7_AQ = POS1_IONS[:(len(POS1_IONS)-2-1)]  # group 1

MOLAR_MASSES = {
    "H": 1.01, "Li": 6.94, "Na": 22.99, "K": 39.10, "Rb": 85.47, "Cs": 132.91, "Fr": 223, "NH4": 17.04,
    "Be": 9.01, "Mg": 24.31, "Ca": 40.08, "Sr": 87.62, "Ba": 137.33, "Ra": 226.00,
    "Sc": 44.96, "Y": 88.91, "La": 138.91, "Ac": 227,
    # lanthanides, actinides to go here
    "Ti": 47.87, "Zr": 91.22, "Hf": 178.49, "Rf": 261,
    "V": 50.94, "Nb": 92.91, "Ta": 180.95,
    "Cr": 52.00, "Mo": 95.94, "W": 183.84,
    "Mn": 54.94, "Tc": 98, "Re": 186.21,
    "Fe": 55.85, "Ru": 101.07, "Os": 190.23,
    "Co": 58.93, "Rh": 102.91, "Ir": 192.22,
    "Ni": 58.69, "Pd": 106.42, "Pt": 195.08,
    "Cu": 63.55, "Ag": 107.87, "Au": 196.97,
    "Zn": 65.41, "Cd": 112.41, "Hg": 200.59, "Hg2": 401.18,
    "B": 10.81, "Al": 26.98, "Ga": 69.72, "In": 114.82, "Tl": 204.38,
    "C": 12.01, "Si": 28.09, "Ge": 72.64, "Sn": 118.71, "Pb": 207.2,
    "N": 14.01, "P": 30.97, "As": 74.92, "Sb": 121.76, "Bi": 208.98,
    "O": 16.00, "S": 32.07, "Se": 78.96, "Te" 127.60, "Po": 209,
    "F": 19.00, "Cl": 35.45, "Br": 79.90, "I": 126.90, "At": 210
}
OUTPUT_PHRASE = [["The mass of", "product", "is"], ["The", "two", "ions", "will", "not", "create", "a", "precipitate."]]


# ---------- FUNCTIONS ---------- #
def startingScreen():
    # wasn't sure what the starting screen was supposed to be, but it was mentioned in the requirements so
    print("""
Welcome! This precipitate calculator can:
- tell whether 2 ions form a precipitate
- calculate the mass of the precipitate if one forms
    """)


# ----- Inputs
def getReactants():
    """
    obtain reactant's element symbol, volume of solution, conc. of solution
    :return: list (--> str)
    """
    ion_info = [[], []]
    charge = {0: "positive", 1: "negative"}
    for i in range(2):
        # --- INPUTS
        ion = input(f"Please enter the {charge[i]} reacting ion. (no charges) ")
        volume = input(f"What is the volume of the {charge[i]} solution? (L) ")
        # check
        conc = input(f"What is the concentration of the {charge[i]} solution? (mol/L) ")
        # check

        # --- PROCESSING
        # adding values
        ion_info[i].append(ion)
        ion_info[i].append(volume)
        ion_info[i].append(conc)
        # modifying values: str to float for volume and concentration values
        for j in range(1, 3):
            str_to_float_value = float(ion_info[i][j])
            ion_info[i][j] = str_to_float_value

    # --- OUTPUT
    return ion_info


# ----- Processing
def determinePrecipitate(cation, anion):
    global COL2_PRECIPITATES, COL3_PRECIPITATES, COL4_PRECIPITATES, COL5_TO_7_AQ
    # column 2
    if (anion == "ClO4" and cation in ("Rb", "Cs")) or (anion == "CH#COO" and cation in ("Ag", "Hg2")):
        return True
    # column 3-7
    elif anion == "F" and cation in COL2_PRECIPITATES:
        return True
    elif anion in ("Cl", "Br", "I") and cation in COL3_PRECIPITATES:
        return True
    elif anion == "SO4" and cation in COL4_PRECIPITATES:
        return True
    elif anion in ("CO3", "PO4", "SO3") and cation not in COL5_TO_7_AQ:
        return True
    elif anion in ("IO3", "OOCCOO") and (cation not in COL5_TO_7_AQ or cation not in ("NH4", "Co")):
        return True
    elif anion == "OH" and cation not in COL5_TO_7_AQ:
        return True
    else:
        return False


def getCharge(elem):
    global POS1_IONS, POS2_IONS, POS3_IONS, POS4_IONS, POS5_IONS, POS6_IONS, POS7_IONS
    global NEG1_IONS, NEG2_IONS, NEG3_IONS, NEG4_IONS
    if elem in POS1_IONS:
        return 1
    elif elem in POS2_IONS:
        return 2
    elif elem in POS3_IONS:
        return 3
    elif elem in POS4_IONS:
        return 4
    elif elem in POS5_IONS:
        return 5
    elif elem in POS6_IONS:
        return 6
    elif elem in POS7_IONS:
        return 7
    elif elem in NEG1_IONS:
        return -1
    elif elem in NEG2_IONS:
        return -2
    elif elem in NEG3_IONS:
        return -3
    else:
        return -4


def getMolarMass(elem):
    global MOLAR_MASSES
    return MOLAR_MASSES[elem]


def GCF(a, b):
    if b == 0:
        return a
    return GCF(b, a % b)


def LCM(a, b):
    g = GCF(a, b)
    return a // g * b


def balanceProduct(cation, anion, cation_charge, anion_charge):
    charge_magnitude = LCM(cation_charge, abs(anion_charge))
    # neither subscript is going to have remainder b/c LCM of both was found (definitely divisible)
    cation_subscript = charge_magnitude // cation_charge
    anion_subscript = charge_magnitude // abs(anion_charge)
    # product = [cation, cation_subscript, anion, anion_subscript]
    if cation_subscript == 1 and anion_subscript == 1:
        product = f"{cation}{anion}"
    elif cation_subscript == 1 and anion_subscript != 1:
        product = f"{cation}{anion}{anion_subscript}"
    elif cation_subscript != 1 and anion_subscript == 1:
        product = f"{cation}{cation_subscript}{anion}"
    else:
        product = f"{cation}{cation_subscript}{anion}{anion_subscript}"

    return product, cation_subscript, anion_subscript


def findNumAtoms(elem):  # only if diatomics need to be considered
    last_char = elem[-1]
    if last_char.isnumeric():
        return int(last_char)
    else:
        return 1


def balanceEquation(product_cation_subscript, product_anion_subscript):
    cation_atoms, anion_atoms = 1, 1
    # product_cation_atoms = product_cation_subscript
    # product_anion_atoms = product_anion_subscript

    min_cation_atoms = LCM(cation_atoms, product_cation_subscript)
    min_anion_atoms = LCM(anion_atoms, product_anion_subscript)

    cation_coeff = min_cation_atoms // cation_atoms
    anion_coeff = min_anion_atoms // anion_atoms

    if product_cation_subscript > product_anion_subscript:
        product_coeff = min_cation_atoms // product_cation_subscript
    else:
        product_coeff = min_anion_atoms // product_anion_subscript

    return cation_coeff, anion_coeff, product_coeff


def findLRMoles(ion_info_set, cation_coeff, anion_coeff):
    cation_volume = ion_info_set[0][1]
    cation_conc = ion_info_set[0][2]
    cation_moles = cation_conc * cation_volume

    anion_volume = ion_info_set[1][1]
    anion_conc = ion_info_set[1][2]
    anion_moles = anion_conc * anion_volume

    # assume cation is LR
    theoretical_anion_moles = cation_moles * anion_coeff / cation_coeff
    if theoretical_anion_moles > anion_moles:  # return moles of LR
        print(f"The limiting reagent is {ion_info_set[1][0]}.")
        return anion_coeff, anion_moles  # anion is LR
    else:
        print(f"The limiting reagent is {ion_info_set[0][0]}.")
        return cation_coeff, cation_moles  # cation is LR


def findPrecipitateMoles(LR_moles, LR_coeff, product_coeff):
    product_moles = LR_moles * product_coeff / LR_coeff
    return product_moles


def calcPrecipitateMass(cation_mass, anion_mass, product_moles, cation_subscript, anion_subscript):
    product_molar_mass = cation_mass * cation_subscript + anion_mass * anion_subscript  # in g/mol
    product_mass = product_molar_mass * product_moles
    return f"{product_mass:.2f}"



# ----- Outputs
def outputPrecipitateCalcResult(precipitate_boolean, product=None, precipitate_mass=None):
    global OUTPUT_PHRASE
    if precipitate_boolean:  # if a precipitate does form
        # modify value in subarray
        OUTPUT_PHRASE[0].pop(1)
        OUTPUT_PHRASE[0].insert(1, product)
        # add value to subarray
        OUTPUT_PHRASE[0].append(precipitate_mass)
        OUTPUT_PHRASE[0].append("grams.")

        print(*OUTPUT_PHRASE[0])
    else:
        print(*OUTPUT_PHRASE[1])


def again():
    again = input("Use the calculator again? (y/N) ")
    if again.lower() == "y" or again.lower() == "yes":
        pass
    else:
        exit()


# ---------- MAIN PROGRAM CODE ---------- #
if __name__ == "__main__":
    while True:
        startingScreen()
        ### INPUTS
        IONS = getReactants()
        CATION, ANION = IONS[0][0], IONS[1][0]

        ### PROCESSING
        PRECIPITATE = determinePrecipitate(CATION, ANION)
        if PRECIPITATE:  # if precipitate forms
            CATION_CHARGE = getCharge(CATION)
            ANION_CHARGE = getCharge(ANION)
            CATION_MASS = getMolarMass(CATION)
            ANION_MASS = getMolarMass(ANION)

            PRODUCT_CHEM_FORMULA, CATION_SUBSCRIPT, ANION_SUBSCRIPT = balanceProduct(CATION, ANION, CATION_CHARGE, ANION_CHARGE)

            # balanced equation to obtain coefficients to use for mole ratio calcs
            CATION_COEFF, ANION_COEFF, PRODUCT_COEFF = balanceEquation(CATION_SUBSCRIPT, ANION_SUBSCRIPT)

            # mole ratio calculations
            LR_COEFF, LR_MOLES = findLRMoles(IONS, CATION_COEFF, ANION_COEFF)
            PRECIPITATE_MOLES = findPrecipitateMoles(LR_MOLES, LR_COEFF, PRODUCT_COEFF)

            # calculate precipitate mass
            PRECIPITATE_MASS = calcPrecipitateMass(CATION_MASS, ANION_MASS, PRECIPITATE_MOLES, CATION_SUBSCRIPT, ANION_SUBSCRIPT)

            ### OUTPUT
            outputPrecipitateCalcResult(PRECIPITATE, PRODUCT_CHEM_FORMULA, PRECIPITATE_MASS)
        else:
            print(*OUTPUT_PHRASE[1])
        again()
