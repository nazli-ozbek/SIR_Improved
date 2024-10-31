import numpy as np
import matplotlib.pyplot as plt

def SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks):
    # Initialize arrays for populations (SIR for A, B, AB, and BA) with np.empty
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

    Sb[0] = Nb * 0.99
    Ib[0] = Nb * 0.01
    Rb_[0] = 0

    # Simulation loop
    for t in range(1, weeks + 1):
        # Visitor inflows and outflows, all groups now included
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

        # Update group A (local population)
        Sa[t] = max(0, min(Na, Sa[t - 1] - Ka * Sa[t - 1] * Ia[t - 1] + inflow_Sba - outflow_Sa))
        Ia[t] = max(0, min(Na - Sa[t], Ia[t - 1] + Ka * Sa[t - 1] * Ia[t - 1] - Ra * Ia[t - 1] + inflow_Iba - outflow_Ia))
        Ra_[t] = max(0, min(Na - Sa[t] - Ia[t], Ra_[t - 1] + Ra * Ia[t - 1] + inflow_Rba - outflow_Ra))


        # Update group B (local population)
        Sb[t] = max(0, min(Nb, Sb[t - 1] - Kb * Sb[t - 1] * Ib[t - 1] + inflow_Sab - outflow_Sb))
        Ib[t] = max(0, min(Nb - Sb[t], Ib[t - 1] + Kb * Sb[t - 1] * Ib[t - 1] - Rb * Ib[t - 1] + inflow_Iab - outflow_Ib))
        Rb_[t] = max(0, min(Nb - Sb[t] - Ib[t], Rb_[t - 1] + Rb * Ib[t - 1] + inflow_Rab - outflow_Rb))


        # Update visitors from A in B (AB)
        Sab[t] = max(0, inflow_Sab - Kb * Sab[t - 1] * Ib[t - 1])
        Iab[t] = max(0, inflow_Iab + Kb * Sab[t - 1] * Ib[t - 1] - Rb * Iab[t - 1])
        Rab[t] = max(0, inflow_Rab + Rb * Iab[t - 1])


        # Update visitors from B in A (BA)
        Sba[t] = max(0, inflow_Sba - Ka * Sba[t - 1] * Ia[t - 1])
        Iba[t] = max(0, inflow_Iba + Ka * Sba[t - 1] * Ia[t - 1] - Ra * Iba[t - 1])
        Rba[t] = max(0, inflow_Rba + Ra * Iba[t - 1])

        print(Sa[t]+Ia[t]+Ra_[t]+ Sb[t]+Ib[t]+Rb_[t]+Sab[t]+Iab[t]+Rab[t]+Sba[t]+Iba[t]+Rba[t])

    return (Sa, Ia, Ra_, Sb, Ib, Rb_, Sab, Iab, Rab, Sba, Iba, Rba)

if __name__ == "__main__":
    # Parameters
    Na, Nb = 1000, 800  # Population sizes of groups A and B
    Ka, Kb = 0.002, 0.001  # Transmission coefficients for A and B
    Ra, Rb = 0.1, 0.1  # Recovery rates for A and B
    Mab, Mba = 0.05, 0.03  # Movement rates between A and B
    Dab, Dba = 2, 3  # Average stay durations in weeks
    weeks = 100  # Duration of the simulation

    # Run the SIR model
    results = SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks)

    # Generate time axis for plotting
    time = np.arange(weeks + 1)

    # Plot the results
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))  # 2x2 grid, larger figure size

    # Plot 1: Group A (Local Population)
    axs[0, 0].plot(time, results[0], label='Susceptible A', color='green', linewidth=2)
    axs[0, 0].plot(time, results[1], label='Infected A', color='red', linewidth=2)
    axs[0, 0].plot(time, results[2], label='Recovered A', color='blue', linewidth=2)
    axs[0, 0].set_title('Group A (Local Population)', fontsize=16)
    axs[0, 0].legend(fontsize=10)
    axs[0, 0].tick_params(axis='both', which='major', labelsize=10)

    # Plot 2: Group B (Local Population)
    axs[0, 1].plot(time, results[3], label='Susceptible B', color='green', linewidth=2)
    axs[0, 1].plot(time, results[4], label='Infected B', color='red', linewidth=2)
    axs[0, 1].plot(time, results[5], label='Recovered B', color='blue', linewidth=2)
    axs[0, 1].set_title('Group B (Local Population)', fontsize=16)
    axs[0, 1].legend(fontsize=10)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=10)

    # Plot 3: Visitors from A in B (AB)
    axs[1, 0].plot(time, results[6], label='Susceptible AB', color='orange', linewidth=2)
    axs[1, 0].plot(time, results[7], label='Infected AB', color='purple', linewidth=2)
    axs[1, 0].plot(time, results[8], label='Recovered AB', color='cyan', linewidth=2)
    axs[1, 0].set_title('Visitors from A in B (AB)', fontsize=16)
    axs[1, 0].legend(fontsize=10)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=10)

    # Plot 4: Visitors from B in A (BA)
    axs[1, 1].plot(time, results[9], label='Susceptible BA', color='orange', linewidth=2)
    axs[1, 1].plot(time, results[10], label='Infected BA', color='purple', linewidth=2)
    axs[1, 1].plot(time, results[11], label='Recovered BA', color='cyan', linewidth=2)
    axs[1, 1].set_title('Visitors from B in A (BA)', fontsize=16)
    axs[1, 1].legend(fontsize=10)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=10)

    # Final adjustments
    for ax in axs.flat:
        ax.set_xlabel('Weeks', fontsize=12)
        ax.set_ylabel('Population', fontsize=12)

    plt.tight_layout()
    plt.show()
