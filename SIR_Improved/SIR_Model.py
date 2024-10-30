import numpy as np
import matplotlib.pyplot as plt

#Na & Nb >> population of two groups
#Sa & Ia & Ra >> SIR for A
#Sb & Ib & Rb >> SIR for B
#Sab & Iab & Rab >> SIR for AB
#Sba & Iba & Rba >> SIR for BA
#Ka & Kb >> transmission coefficients for group A and B (will be taken different)
#Ra & Rb >> recovery rates for A and B
#Mab >> the proportion of A that travels to B / stay during Dab weeks
#Mba >> the proportion of B that travels to A / stay during Dba weeks
#weeks >> duration of the disease

def SIR_Model(Na, Nb, Ia0, Ib0, Ra, Rb, Ka, Kb, Mab, Mba, weeks):
    #initial SIR for A
    Sa = np.empty(weeks + 1)
    Ia = np.empty(weeks + 1)
    Ra = np.empty(weeks + 1)

    Sa[0] = Na - Ia0
    Ia[0] = Ia0
    Ra[0] = 0

    #initial SIR for B
    Sb = np.empty(weeks + 1)
    Ib = np.empty(weeks + 1)
    Rb = np.empty(weeks + 1)

    Sb[0] = Nb - Ib0
    Ib[0] = Ib0
    Rb[0] = 0

    #initial SIR for A in B
    Sab = np.empty(weeks + 1)
    Iab = np.empty(weeks + 1)
    Rab = np.empty(weeks + 1)
    Sab[0] = 0
    Iab[0] = 0
    Rab[0] = 0

    #initial SIR for B in A
    Sba = np.empty(weeks + 1)
    Iba = np.empty(weeks + 1)
    Rba = np.empty(weeks + 1)
    Sba[0] = 0
    Iba[0] = 0
    Rba[0] = 0

    # for t in range(weeks):





if __name__ == "__main__":
    SIR_Model(5, 1000, 24, 5/3, 0.00140704)