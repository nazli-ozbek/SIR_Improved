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

def SIR_Model(Na, Nb, Ia0, Ib0, Reca, Recb, Ka, Kb, Mab, Mba, weeks, d):
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

    for t in range(weeks):
        #SIR for A
        Sa = Na - Ia[t] - Ra[t] - Sab[t-1]
        Ia = Ia[t-1] - Reca * Ia[t-1] + Ka * Sa[t-1] * (Ia[t-1] + Iab[t-1])
        Ra = Ra[t-1] + Reca * Ia[t-1]

        # SIR for B
        Sb = Nb - Ib[t] - Rb[t] - Sba[t - 1]
        Ib = Ib[t - 1] - Recb * Ib[t - 1] + Kb * Sb[t - 1] * (Ib[t - 1] + Iba[t - 1])
        Rb = Rb[t - 1] + Recb * Ib[t - 1]

        #SIR for A to B
        Sab[t] = Sab[t-1] + Sa[t] * Mab[t]
        Iab[t] = Iab[t-1] - Rab[t-1] + Ia[t] * Mab
        Rab[t] = Rab[t-1] + Ra[t] * Mab

        # SIR for B to A
        Sba[t] = Sba[t - 1] + Sb[t] * Mba[t]
        Iba[t] = Iba[t - 1] - Rba[t - 1] + Ib[t] * Mba
        Rba[t] = Rba[t - 1] + Rb[t] * Mba







if __name__ == "__main__":
    SIR_Model(5, 1000, 24, 5/3, 0.00140704)