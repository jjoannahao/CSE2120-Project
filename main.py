"""
title: Precipitate Calculator
author: Joanna Hao
date-created: 2022-10-20
"""
# ---------- VARIABLES ---------- #
POS1_IONS = ("H", "Li", "Na", "K", "Rb", "Cs", "Fr", "NH4", "Ag", "Tl")
POS2_IONS = ("Be", "Mg", "Ca", "Sr", "Ba", "Ra", "Mn", "Co", "Ni", "Pd", "Cu", "Zn", "Cd", "Pb", "Po")

MOLAR_MASSES = {
    "H": 1.01, "Li": 6.94, "Na": 22.99, "K": 39.10, "Rb": 85.47, "Cs": 132.91, "Fr": 223,


    "Pb": 207.2,

    "F": 19.00, "Cl": 35.45, "Br": 79.90, "I": 126.90, "At": 210
}

# ---------- FUNCTIONS ---------- #
# ----- Inputs
def getReactants():
    """
    obtain reactant's element symbol, volume of solution, conc. of solution
    :return: list (--> str, int)
    """
    ion_info = []
    charge = {0: "positive", 1: "negative"}
    for i in range(2):
        ion = input(f"Please enter the {charge[i]} reacting ion. (no charges) ")
        volume = input(f"What is the volume of the {charge[i]} solution? (L) ")
        # check
        conc = input(f"What is the concentration of the {charge[i]} solution? (mol/L) ")
        # check
        ion_info.append((ion, float(volume), float(conc)))
    return ion_info


# ----- Processing
def determinePrecipitate(cation, anion):
    # column 1

    # column 2
    if anion == "F" and cation in ("Li", "Mg", "Ca", "Sr", "Ba", "Fe", "Hg2", "Pb"):
        return True
    elif anion in ("Cl", "Br", "I") and cation in ("Ag", "Hg2", "Pb", "Tl", "Cu"): # Cu has to be less common +1 charge
        return True
    elif



def getCharge(elem):
    global POS1_IONS, POS2_IONS  # , all the remaining charges pos & neg
    if elem in POS1_IONS:
        return 1
    elif elem in POS2_IONS:
        return 2


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
    """
    determine balanced chemical formula of product
    :param cation: str
    :param anion: str
    :param cation_charge: int
    :param anion_charge: int
    :return:
    """
    # FINDING PRODUCT
    charge_magnitude = LCM(cation_charge, abs(anion_charge))
    # neither subscript is going to have remainder b/c LCM w/ other was found
    cation_subscript = charge_magnitude // cation_charge
    anion_subscript = charge_magnitude // abs(anion_charge)
    product = f"{cation}{cation_subscript}{anion}{anion_subscript}"

    # FINDING COEFFICIENTS
    # check if anion is diatomic (BrINClHOF)/polyatomic (P4, S8)
    diatomics = ("Br", "I", "N", "Cl", "H", "O", "F")
    polyatomics = {"P": "P4", "S": "S8"}
    if anion in diatomics:
        anion = f"{anion}2"
    elif anion in polyatomics:
        anion = polyatomics[anion]

    coeff1 = cation_subscript  # cation always a single atom -- matches whatever the product cation subscript it
    coeff2 = LCM(int(anion[-1]), anion_subscript) // int(anion[-1])
    coeff3 = LCM(int(anion[-1]), anion_subscript) // anion_subscript
    chem_equation = [coeff1, cation, coeff2, anion, coeff3, product]
    return chem_equation


def findLR(volume, conc, elem1_coeff, elem2_coeff):
    # use LR cases
    pass


def findMoles():
    # use LR, mole ratio
    pass


def calcMass():
    # find mass of product (use molar mass, moles)
    # only if a precipitate is formed, otherwise print smth like 'doesn't form a precipitate'
    pass


# ----- Outputs



# ---------- MAIN PROGRAM CODE ---------- #
if __name__ == "__main__":
    IONS = getReactants()
    # print(IONS)
    # print(getCharge("OH"))
