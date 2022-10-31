# CSE2120-Project

Precipitate Calculator

# Info from Document
## Data Structures, Functions Required
* tuple
* 1D array 
  * tuple crosses this off or might need another one
* dict optional
- +1 functions must 
  - add value to array
  - modify a value in an array
  - final output (list) -- add mass value, modify placeholder of product chemical formula 

### DONE
* 2D array (list/tuple inside list/tuple)
  * 3 items in broader list -- 2 reactants, 1 product
    * each sublist has all necessary properties of substances
      * could be substances = ((cation, vol, conc), (anion, vol, conc), (product))


## Document's Pseudocode
inputs include:
* name of ion, volume, conc.
   * need to know charge, molar mass

processing
* determine which product is the precipitate (all formation rxns?)
* balance equation 
  * use a check to make sure net charge (of all things) is 0
* find LR
  * use volume, conc., coefficients from chem rxn
* find moles of product/precipitate
  * use LR, mole ratio (coefficients)
* find (theoretical) molar mass of precipitate
* calculate mass of precipitate w/ molar mass, moles

outputs
* mass of precipitate