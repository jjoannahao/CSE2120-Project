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
             "Mn", "Co", "Ni", "Pd", "Cu", "Zn", "Cd", "Pb", "Po")
NEG1_IONS = ("F", "Cl", "Br", "I", "At",
             "CH3COO", "C6H5COO", "HCO3", "ClO4", "ClO3", "ClO2", "OCl", "ClO", "CN", "OH", "IO3", "NO3", "NO2", "HOOCCOO", "MnO4", "H2PO4", "HSO4", "HSO3", "HS", "SCN")
NEG2_IONS = ("O", "S", "Se", "Te",
             "C2", "CO3", "CrO4", "Cr2O7", "OOCCOO", "O2", "S2", "HPO4", "SiO3", "SO4", "SO3", "S2O3")

COL2_PRECIPITATES = ("Li", "Mg", "Ca", "Sr", "Ba", "Fe", "Hg2", "Pb")
COL3_PRECIPITATES = ("Cu", "Ag", "Hg2", "Pb", "Tl")  # Cu = +1 charge (less common)
COL4_PRECIPITATES = ("Ca", "Sr", "Ba", "Ag", "Hg2", "Pb", "Ra")
COL5_TO_7_AQ = POS1_IONS[:(len(POS1_IONS)-2-1)]  # group 1


# not in solubility table
# NEG3_IONS = ("N", "P", "As",
#              "BO3", "PO4", )

MOLAR_MASSES = {
    "H": 1.01, "Li": 6.94, "Na": 22.99, "K": 39.10, "Rb": 85.47, "Cs": 132.91, "Fr": 223, "NH4": 17.04,


    "Pb": 207.2,

    "F": 19.00, "Cl": 35.45, "Br": 79.90, "I": 126.90, "At": 210
}

# ---------- FUNCTIONS ---------- #
def startingScreen():
    pass


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

    if (anion == "ClO4" and cation in ("Rb", "Cs")) or (anion == "CH#COO" and cation in ("Ag", "Hg2")):  # column 2
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
    global POS1_IONS, POS2_IONS, NEG1_IONS, NEG2_IONS
    if elem in POS1_IONS:
        return 1
    elif elem in POS2_IONS:
        return 2
    elif elem in NEG1_IONS:
        return -1
    else:
        return -2


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
    product = [cation, cation_subscript, anion, anion_subscript]
    # product = f"{cation}{cation_subscript}{anion}{anion_subscript}"
    return product


def findNumAtoms(elem):  # only if diatomics need to be considered
    last_char = elem[-1]
    if last_char.isnumeric():
        return int(last_char)
    else:
        return 1


def balanceEquation(product):
    cation_atoms, anion_atoms = 1, 1
    product_cation_atoms = product[1]
    product_anion_atoms = product[3]

    min_cation_atoms = LCM(cation_atoms, product_cation_atoms)
    min_anion_atoms = LCM(anion_atoms, product_anion_atoms)

    cation_coeff = min_cation_atoms // cation_atoms
    anion_coeff = min_anion_atoms // anion_atoms

    if product_cation_atoms > product_anion_atoms:
        product_coeff = min_cation_atoms // product_cation_atoms
    else:
        product_coeff = min_anion_atoms // product_anion_atoms

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
        return anion_coeff, anion_moles  # anion is LR
    else:
        return cation_coeff, cation_moles  # cation is LR


def findPrecipitateMoles(LR_moles, LR_coeff, product_coeff):
    product_moles = LR_moles * product_coeff / LR_coeff
    return product_moles


def calcPrecipitateMass(cation_mass, anion_mass, product_moles, cation_coeff, anion_coeff):
    product_molar_mass = cation_mass * cation_coeff + anion_mass * anion_coeff  # in g/mol
    product_mass = product_molar_mass * product_moles
    return product_mass



# ----- Outputs
def outputPrecipitateCalcResult(precipitate_boolean, product=None, precipitate_mass=None):
    output_phrase = [["The mass of", "product", "is"], ["The", "two", "ions", "will", "not", "form", "a", "precipitate."]]

    if precipitate_boolean:  # if a precipitate does form
        # modify value in subarray
        output_phrase[0].pop(1)
        output_phrase[0].insert(1, product)
        # add value to subarray
        output_phrase[0].append(precipitate_mass)
        output_phrase[0].append("grams.")

        print(*output_phrase[0])
    else:
        print(*output_phrase[1])


# ---------- MAIN PROGRAM CODE ---------- #
if __name__ == "__main__":
    ### INPUTS
    IONS = getReactants()

    ### PROCESSING
    CATION, ANION = IONS[0][0], IONS[1][0]
    PRECIPITATE = determinePrecipitate(CATION, ANION)
    if PRECIPITATE:  # if precipitate forms
        # obtaining charges, molar masses of reactants
        CATION_CHARGE = getCharge(CATION)
        ANION_CHARGE = getCharge(ANION)
        CATION_MASS = getMolarMass(CATION)
        ANION_MASS = getMolarMass(ANION)

        # obtained product (list of cation, subsript, anion, subscript)
        PRODUCT = balanceProduct(CATION, ANION, CATION_CHARGE, ANION_CHARGE)
        PRODUCT_CHEM_FORMULA = " ".join(PRODUCT)
        # balanced equation to obtain coefficients to use for mole ratio calcs
        CATION_COEFF, ANION_COEFF, PRODUCT_COEFF = balanceEquation(PRODUCT)

        # mole ratio calculations
        LR_COEFF, LR_MOLES = findLRMoles(IONS, CATION_COEFF, ANION_COEFF)
        PRECIPITATE_MOLES = findPrecipitateMoles(LR_MOLES, LR_COEFF, PRODUCT_COEFF)

        # calculate precipitate mass
        PRECIPITATE_MASS = calcPrecipitateMass(CATION_MASS, ANION_MASS, PRECIPITATE_MOLES, CATION_COEFF, ANION_COEFF)

    ### OUTPUT
    outputPrecipitateCalcResult(PRECIPITATE, PRODUCT_CHEM_FORMULA, PRECIPITATE_MASS)
