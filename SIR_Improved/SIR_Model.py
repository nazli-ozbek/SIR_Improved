import numpy as np
import matplotlib.pyplot as plt

# Na & Nb >> population of two groups
# Sa & Ia & Ra >> SIR for A
# Sb & Ib & Rb >> SIR for B
# Sab & Iab & Rab >> SIR for AB
# Sba & Iba & Rba >> SIR for BA
# Ka & Kb >> transmission coefficients for group A and B (will be taken different)
# Reca & Recb >> recovery rates for A and B
# Mab >> the proportion of A that travels to B / stay during Dab weeks
# Mba >> the proportion of B that travels to A / stay during Dba weeks
# weeks >> duration of the disease

def SIR_Model(Na, Nb, Ka, Kb, Reca, Recb, Mab, Mba, Dab, Dba, weeks):
    #Creation of arrays
    Sa = np.empty(weeks + 1)
    Ia = np.empty(weeks + 1)
    Ra = np.empty(weeks + 1)

    Sb = np.empty(weeks + 1)
    Ib = np.empty(weeks + 1)
    Rb = np.empty(weeks + 1)

    Sab = np.empty(weeks + 1)
    Iab = np.empty(weeks + 1)
    Rab = np.empty(weeks + 1)

    Sba = np.empty(weeks + 1)
    Iba = np.empty(weeks + 1)
    Rba = np.empty(weeks + 1)

    #Initialization of SIR of A and B
    Sa[0] = Na * 0.99
    Ia[0] = Na * 0.01
    Ra[0] = 0

    Sb[0] = Nb
    Ib[0] = 0
    Rb[0] = 0

    for t in range(1, weeks + 1):
        #Calculation of inflow and outflows
        inflow_Sab = Mab * Sa[t - 1] / Dab
        inflow_Iab = Mab * Ia[t - 1] / Dab
        inflow_Rab = Mab * Ra[t - 1] / Dab

        inflow_Sba = Mba * Sb[t - 1] / Dba
        inflow_Iba = Mba * Ib[t - 1] / Dba
        inflow_Rba = Mba * Rb[t - 1] / Dba

        outflow_Sa = Mab * Sa[t - 1] / Dab
        outflow_Ia = Mab * Ia[t - 1] / Dab
        outflow_Ra = Mab * Ra[t - 1] / Dab

        outflow_Sb = Mba * Sb[t - 1] / Dba
        outflow_Ib = Mba * Ib[t - 1] / Dba
        outflow_Rb = Mba * Rb[t - 1] / Dba

        #Application of difference equations
        # A
        Sa[t] = max(0, min(Na, Sa[t - 1] - Ka * Sa[t - 1] * Ia[t - 1] + inflow_Sba - outflow_Sa))
        Ia[t] = max(0, min(Na - Sa[t] - Ra[t], Ia[t - 1] + Ka * Sa[t - 1] * Ia[t - 1] - Reca * Ia[t - 1] + inflow_Iba - outflow_Ia))
        Ra[t] = max(0, min(Na - Sa[t] - Ia[t], Ra[t - 1] + Reca * Ia[t - 1] + inflow_Rba - outflow_Ra))

        # B
        Sb[t] = max(0, min(Nb, Sb[t - 1] - Kb * Sb[t - 1] * Ib[t - 1] + inflow_Sab - outflow_Sb))
        Ib[t] = max(0, min(Nb - Sb[t] - Rb[t], Ib[t - 1] + Kb * Sb[t - 1] * Ib[t - 1] - Recb * Ib[t - 1] + inflow_Iab - outflow_Ib))
        Rb[t] = max(0, min(Nb - Sb[t] - Ib[t], Rb[t - 1] + Recb * Ib[t - 1] + inflow_Rab - outflow_Rb))

        # AB
        Sab[t] = max(0, inflow_Sab - Kb * Sab[t - 1] * Ib[t - 1])
        Iab[t] = max(0, inflow_Iab + Kb * Sab[t - 1] * Ib[t - 1] - Recb * Iab[t - 1])
        Rab[t] = max(0, inflow_Rab + Recb * Iab[t - 1])

        # BA
        Sba[t] = max(0, inflow_Sba - Ka * Sba[t - 1] * Ia[t - 1])
        Iba[t] = max(0, inflow_Iba + Ka * Sba[t - 1] * Ia[t - 1] - Reca * Iba[t - 1])
        Rba[t] = max(0, inflow_Rba + Reca * Iba[t - 1])


    #Plotting
    time = np.arange(weeks + 1)

    fig, axs = plt.subplots(2, 2, figsize=(18, 12))

    # Plot 1: Group A (Local Population)
    axs[0, 0].plot(time, Sa, label='Susceptible', color='blue')
    axs[0, 0].plot(time, Ia, label='Infected', color='red')
    axs[0, 0].plot(time, Ra, label='Recovered', color='green')
    axs[0, 0].set_title('Group A')
    axs[0, 0].set_xlabel('Weeks')
    axs[0, 0].set_ylabel('Population')
    axs[0, 0].legend()

    # Plot 2: Group B (Local Population)
    axs[0, 1].plot(time, Sb, label='Susceptible', color='blue')
    axs[0, 1].plot(time, Ib, label='Infected', color='red')
    axs[0, 1].plot(time, Rb, label='Recovered', color='green')
    axs[0, 1].set_title('Group B')
    axs[0, 1].set_xlabel('Weeks')
    axs[0, 1].set_ylabel('Population')
    axs[0, 1].legend()

    # Plot 3: Visitors from A in B
    axs[1, 0].plot(time, Sab, label='Susceptible', color='blue')
    axs[1, 0].plot(time, Iab, label='Infected', color='red')
    axs[1, 0].plot(time, Rab, label='Recovered', color='green')
    axs[1, 0].set_title('A in B')
    axs[1, 0].set_xlabel('Weeks')
    axs[1, 0].set_ylabel('Population')
    axs[1, 0].legend()

    # Plot 4: Visitors from B in A
    axs[1, 1].plot(time, Sba, label='Susceptible', color='blue')
    axs[1, 1].plot(time, Iba, label='Infected', color='red')
    axs[1, 1].plot(time, Rba, label='Recovered', color='green')
    axs[1, 1].set_title('B in A')
    axs[1, 1].set_xlabel('Weeks')
    axs[1, 1].set_ylabel('Population')
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    Na, Nb = 1000, 800
    Ka, Kb = 0.002, 0.001
    Reca, Recb = 0.1, 0.05
    Mab, Mba = 0.05, 0.05
    Dab, Dba = 3, 3
    weeks = 16

    SIR_Model(Na, Nb, Ka, Kb, Reca, Recb, Mab, Mba, Dab, Dba, weeks)
