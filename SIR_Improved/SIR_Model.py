import numpy as np
import matplotlib.pyplot as plt

def SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks):
    # Initialize populations (SIR for A, B, AB, and BA)
    Sa, Ia, Ra_ = [Na * 0.99], [Na * 0.01], [0]   # S, I, R for group A
    Sb, Ib, Rb_ = [Nb], [0], [0]       # S, I, R for group B

    # Visitors from A in B (AB) and from B in A (BA)
    Sab, Iab, Rab = [0], [0], [0]  # A's population visiting B
    Sba, Iba, Rba = [0], [0], [0]  # B's population visiting A

    # Simulation loop
    for n in range(weeks):
        # Calculate visitor inflows and outflows
        inflow_Sab = Mab * Sa[n] / Dab
        inflow_Iab = Mab * Ia[n] / Dab
        inflow_Rab = Mab * Ra_[n] / Dab

        inflow_Sba = Mba * Sb[n] / Dba
        inflow_Iba = Mba * Ib[n] / Dba
        inflow_Rba = Mba * Rb_[n] / Dba

        outflow_Sa = Mab * Sa[n] / Dab
        outflow_Sb = Mba * Sb[n] / Dba

        # Update group A (local population)
        Sa_next = max(0, min(Na, Sa[n] - Ka * Sa[n] * Ia[n] + inflow_Sba - outflow_Sa))
        Ia_next = max(0, min(Na - Sa_next, Ia[n] + Ka * Sa[n] * Ia[n] - Ra * Ia[n]))
        Ra_next = max(0, min(Na - Sa_next - Ia_next, Ra_[n] + Ra * Ia[n] + inflow_Rba))

        # Update group B (local population)
        Sb_next = max(0, min(Nb, Sb[n] - Kb * Sb[n] * Ib[n] + inflow_Sab - outflow_Sb))
        Ib_next = max(0, min(Nb - Sb_next, Ib[n] + Kb * Sb[n] * Ib[n] - Rb * Ib[n]))
        Rb_next = max(0, min(Nb - Sb_next - Ib_next, Rb_[n] + Rb * Ib[n] + inflow_Rab))

        # Update visitors from A in B (AB)
        Sab_next = max(0, inflow_Sab - Kb * Sab[n] * Ib[n])
        Iab_next = max(0, inflow_Iab + Kb * Sab[n] * Ib[n] - Rb * Iab[n])
        Rab_next = max(0, inflow_Rab + Rb * Iab[n])

        # Update visitors from B in A (BA)
        Sba_next = max(0, inflow_Sba - Ka * Sba[n] * Ia[n])
        Iba_next = max(0, inflow_Iba + Ka * Sba[n] * Ia[n] - Ra * Iba[n])
        Rba_next = max(0, inflow_Rba + Ra * Iba[n])

        # Append results for the next step
        Sa.append(Sa_next)
        Ia.append(Ia_next)
        Ra_.append(Ra_next)

        Sb.append(Sb_next)
        Ib.append(Ib_next)
        Rb_.append(Rb_next)

        Sab.append(Sab_next)
        Iab.append(Iab_next)
        Rab.append(Rab_next)

        Sba.append(Sba_next)
        Iba.append(Iba_next)
        Rba.append(Rba_next)

    return (Sa, Ia, Ra_, Sb, Ib, Rb_, Sab, Iab, Rab, Sba, Iba, Rba)

def plot_results(time, Sa, Ia, Ra_, Sb, Ib, Rb_, Sab, Iab, Rab, Sba, Iba, Rba):
    fig, axs = plt.subplots(4, 1, figsize=(12, 20))  # 4 rows, 1 column

    # Plot 1: Group A (Local Population)
    axs[0].plot(time, Sa, label='Susceptible A', color='green')
    axs[0].plot(time, Ia, label='Infected A', color='red')
    axs[0].plot(time, Ra_, label='Recovered A', color='blue')
    axs[0].set_title('Group A (Local Population)')
    axs[0].legend()

    # Plot 2: Group B (Local Population)
    axs[1].plot(time, Sb, label='Susceptible B', color='green')
    axs[1].plot(time, Ib, label='Infected B', color='red')
    axs[1].plot(time, Rb_, label='Recovered B', color='blue')
    axs[1].set_title('Group B (Local Population)')
    axs[1].legend()

    # Plot 3: Visitors from A in B (AB)
    axs[2].plot(time, Sab, label='Susceptible AB', color='orange')
    axs[2].plot(time, Iab, label='Infected AB', color='purple')
    axs[2].plot(time, Rab, label='Recovered AB', color='cyan')
    axs[2].set_title('Visitors from A in B (AB)')
    axs[2].legend()

    # Plot 4: Visitors from B in A (BA)
    axs[3].plot(time, Sba, label='Susceptible BA', color='orange')
    axs[3].plot(time, Iba, label='Infected BA', color='purple')
    axs[3].plot(time, Rba, label='Recovered BA', color='cyan')
    axs[3].set_title('Visitors from B in A (BA)')
    axs[3].legend()

    # Final adjustments
    for ax in axs:
        ax.set_xlabel('Weeks')
        ax.set_ylabel('Number of Individuals')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Parameters
    Na, Nb = 1000, 800             # Population sizes of groups A and B
    Ka, Kb = 0.002, 0.001          # Transmission coefficients for A and B
    Ra, Rb = 0.1, 0.1              # Recovery rates for A and B
    Mab, Mba = 0.05, 0.03          # Movement rates between A and B
    Dab, Dba = 2, 3                # Average stay durations in weeks
    weeks = 500                     # Duration of the simulation

    # Run the SIR model
    results = SIR_Model(Na, Nb, Ka, Kb, Ra, Rb, Mab, Mba, Dab, Dba, weeks)

    # Generate time axis for plotting
    time = np.arange(weeks + 1)

    # Plot the results
    plot_results(time, *results)
