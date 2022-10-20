"""
title: Precipitate Calculator
author: Joanna Hao
date-created: 2022-10-20
"""
# ---------- VARIABLES ---------- #
COL1_IONS = ("H", "Li", "Na", "K", "Rb", "Cs", "Fr", "NH4",   "NO3", "ClO3", "ClO4", "CH3COO")
COL2_ION = "F"
COL3_IONS = ("Cl", "Br", "I")
COL4_ION = "SO4"
COL5_7_IONS = ("CO3", "SO3", "PO4", "OH")
COL6_IONS = ("IO3", "OOCCOO")


# ---------- FUNCTIONS ---------- #
def startScreen():
    """
    display start screen text
    :return:
    """


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
        conc = input(f"What is the concentration of the {charge[i]} solution? (mol/L) ")
        ion_info.append((ion, float(volume), float(conc)))
    return ion_info


# ----- Processing
def getCharge(elem):
    """
    get charge of element
    :param elem: str (from tuple)
    :return:
    """
    global COL1_IONS, COL2_ION, COL3_IONS, COL4_ION, COL5_7_IONS, COL6_IONS

    # problematic: CO3, SO3
    if elem in COL1_IONS[:8]:
        return 1
    elif elem in COL1_IONS[8:]:
        return -1
    elif elem in COL2_ION or elem in COL3_IONS or elem == COL5_7_IONS[-1]:
        return -1
    elif elem in COL4_ION or elem == COL5_7_IONS[:2] or elem == COL6_IONS[-1]:
        return -2
    elif elem == COL5_7_IONS[2]:
        return -3
    else:  # elem == COL6_IONS[0]  which is IO3 --> 1-
        return -1


def balanceEquation(elem1, elem2):
    """
    balance equation so that net charge of reactants, products is 0
    :param ELEM1:
    :param ELEM2:
    :return:
    """




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
    # IONS = getReactants()
    # print(IONS)
    print(getCharge("OH"))
