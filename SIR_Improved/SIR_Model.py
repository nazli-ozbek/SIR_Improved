import numpy as np
import matplotlib.pyplot as plt

# Na & Nb >> population of two groups
# Sa & Ia & Ra >> SIR for A
# Sb & Ib & Rb >> SIR for B
# Sab & Iab & Rab >> SIR for AB
# Sba & Iba & Rba >> SIR for BA
# Ka & Kb >> transmission coefficients for group A and B (will be taken different)
# Ra & Rb >> recovery rates for A and B
# Mab >> the proportion of A that travels to B / stay during Dab weeks
# Mba >> the proportion of B that travels to A / stay during Dba weeks
# weeks >> duration of the disease

def SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks):
    Sa = np.empty(weeks + 1)
    Ia = np.empty(weeks + 1)
    Ra_ = np.empty(weeks + 1)

    Sb = np.empty(weeks + 1)
    Ib = np.empty(weeks + 1)
    Rb_ = np.empty(weeks + 1)

    Sab = np.empty(weeks + 1)
    Iab = np.empty(weeks + 1)
    Rab = np.empty(weeks + 1)

    Sba = np.empty(weeks + 1)
    Iba = np.empty(weeks + 1)
    Rba = np.empty(weeks + 1)

    # Initial conditions
    Sa[0] = Na * 0.99
    Ia[0] = Na * 0.01
    Ra_[0] = 0

    Sb[0] = Nb
    Ib[0] = 0
    Rb_[0] = 0

    # Simulation loop
    for t in range(1, weeks + 1):
        inflow_Sab = Mab * Sa[t - 1] / Dab
        inflow_Iab = Mab * Ia[t - 1] / Dab
        inflow_Rab = Mab * Ra_[t - 1] / Dab

        inflow_Sba = Mba * Sb[t - 1] / Dba
        inflow_Iba = Mba * Ib[t - 1] / Dba
        inflow_Rba = Mba * Rb_[t - 1] / Dba

        outflow_Sa = Mab * Sa[t - 1] / Dab
        outflow_Ia = Mab * Ia[t - 1] / Dab
        outflow_Ra = Mab * Ra_[t - 1] / Dab

        outflow_Sb = Mba * Sb[t - 1] / Dba
        outflow_Ib = Mba * Ib[t - 1] / Dba
        outflow_Rb = Mba * Rb_[t - 1] / Dba

        # A
        Sa[t] = max(0, min(Na, Sa[t - 1] - Ka * Sa[t - 1] * Ia[t - 1] + inflow_Sba - outflow_Sa))
        Ia[t] = max(0, min(Na - Sa[t] - Ra_[t], Ia[t - 1] + Ka * Sa[t - 1] * Ia[t - 1] - Ra * Ia[t - 1] + inflow_Iba - outflow_Ia))
        Ra_[t] = max(0, min(Na - Sa[t] - Ia[t], Ra_[t - 1] + Ra * Ia[t - 1] + inflow_Rba - outflow_Ra))

        # B
        Sb[t] = max(0, min(Nb, Sb[t - 1] - Kb * Sb[t - 1] * Ib[t - 1] + inflow_Sab - outflow_Sb))
        Ib[t] = max(0, min(Nb - Sb[t] - Rb_[t], Ib[t - 1] + Kb * Sb[t - 1] * Ib[t - 1] - Rb * Ib[t - 1] + inflow_Iab - outflow_Ib))
        Rb_[t] = max(0, min(Nb - Sb[t] - Ib[t], Rb_[t - 1] + Rb * Ib[t - 1] + inflow_Rab - outflow_Rb))

        # AB
        Sab[t] = max(0, inflow_Sab - Kb * Sab[t - 1] * Ib[t - 1])
        Iab[t] = max(0, inflow_Iab + Kb * Sab[t - 1] * Ib[t - 1] - Rb * Iab[t - 1])
        Rab[t] = max(0, inflow_Rab + Rb * Iab[t - 1])

        # BA
        Sba[t] = max(0, inflow_Sba - Ka * Sba[t - 1] * Ia[t - 1])
        Iba[t] = max(0, inflow_Iba + Ka * Sba[t - 1] * Ia[t - 1] - Ra * Iba[t - 1])
        Rba[t] = max(0, inflow_Rba + Ra * Iba[t - 1])

        # Print total populations to check for issues
        print(Sa[t] + Ia[t] + Ra_[t] + Sb[t] + Ib[t] + Rb_[t] + Sab[t] + Iab[t] + Rab[t] + Sba[t] + Iba[t] + Rba[t])

    # Generate time axis for plotting
    time = np.arange(weeks + 1)

    # Create plots
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))  # 2x2 grid, larger figure size

    # Plot 1: Group A (Local Population)
    axs[0, 0].plot(time, Sa, label='Susceptible', color='blue')
    axs[0, 0].plot(time, Ia, label='Infected', color='red')
    axs[0, 0].plot(time, Ra_, label='Recovered', color='green')
    axs[0, 0].set_title('Group A')
    axs[0, 0].set_xlabel('Weeks')
    axs[0, 0].set_ylabel('Population')
    axs[0, 0].legend()

    # Plot 2: Group B (Local Population)
    axs[0, 1].plot(time, Sb, label='Susceptible', color='blue')
    axs[0, 1].plot(time, Ib, label='Infected', color='red')
    axs[0, 1].plot(time, Rb_, label='Recovered', color='green')
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

    plt.show()

if __name__ == "__main__":
    Na, Nb = 1000, 800
    Ka, Kb = 0.002, 0.001
    Ra, Rb = 0.1, 0.05
    Mab, Mba = 0.05, 0.05
    Dab, Dba = 3, 3
    weeks = 100

    SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks)
