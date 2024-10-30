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

def SIR_Model(Na, Nb, Ia0, Ib0, Reca, Recb, Ka, Kb, Mab, Mba, weeks):
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

    for t in range(1,weeks + 1):
        #SIR for A
        Ra[t] = Ra[t - 1] + Reca * Ia[t - 1] - Mab * Ra[t -1]
        Ia[t] = Ia[t-1] - Mab*Ia[t-1] + Ka*Sa[t-1]*(Ia[t-1] + Iab[t-1]) - Reca * Ia[t-1]
        Sa[t] = Sa[t-1] - Mab*Sa[t-1] - Ia[t] - Ra[t]
        print("A " , Ra[t] + Ia[t] + Sa[t])

        # SIR for B
        Rb[t] = Rb[t - 1] + Recb * Ib[t - 1] - Mba * Rb[t - 1]
        Ib[t] = Ib[t - 1] - Mba * Ib[t - 1] + Kb * Sb[t - 1] * (Ib[t - 1] + Iba[t - 1]) - Recb * Ib[t - 1]
        Sb[t] = Sb[t - 1] - Mba * Sb[t - 1] - Ib[t] - Rb[t]
        print("B " , Rb[t] + Ib[t] + Sb[t])

        #SIR for A to B
        Rab[t] = Mab * Ra[t] + Rab[t-1] + Iab[t-1] * Reca
        Iab[t] = Mab * Ia[t] + Iab[t-1] - Iab[t-1] * Reca
        Sab[t] =  Sa[t] * Mab + Sab[t - 1]
        print("AB " , Rab[t] + Iab[t] + Sab[t])


        # SIR for B to A
        Rba[t] = Mba * Rb[t] + Rba[t - 1] + Iba[t - 1] * Recb
        Iba[t] = Mba * Ib[t] + Iba[t - 1] - Iba[t - 1] * Recb
        Sba[t] = Sb[t] * Mba + Sba[t - 1]
        print("BA " ,Rba[t] + Iba[t] + Sba[t])

        print("Total", Ra[t] + Ia[t] + Sa[t] + Rb[t] + Ib[t] + Sb[t] + Rab[t] + Iab[t] + Sab[t] + Rba[t] + Iba[t] + Sba[t] )


        # Plotting the results
        plt.figure(figsize=(15, 12))

        # Group A SIR
        plt.subplot(4, 1, 1)
        plt.plot(range(weeks + 1), Sa, 'b-', label='Susceptible A')
        plt.plot(range(weeks + 1), Ia, 'r-', label='Infected A')
        plt.plot(range(weeks + 1), Ra, 'g-', label='Recovered A')
        plt.title('SIR Model for Group A')
        plt.xlabel('Weeks')
        plt.ylabel('Population')
        plt.legend()
        plt.grid(True)

        # Group B SIR
        plt.subplot(4, 1, 2)
        plt.plot(range(weeks + 1), Sb, 'b-', label='Susceptible B')
        plt.plot(range(weeks + 1), Ib, 'r-', label='Infected B')
        plt.plot(range(weeks + 1), Rb, 'g-', label='Recovered B')
        plt.title('SIR Model for Group B')
        plt.xlabel('Weeks')
        plt.ylabel('Population')
        plt.legend()
        plt.grid(True)

        # Travelers from A to B
        plt.subplot(4, 1, 3)
        plt.plot(range(weeks + 1), Sab, 'b-', label='Susceptible in B from A')
        plt.plot(range(weeks + 1), Iab, 'r-', label='Infected in B from A')
        plt.plot(range(weeks + 1), Rab, 'g-', label='Recovered in B from A')
        plt.title('Travelers from A to B')
        plt.xlabel('Weeks')
        plt.ylabel('Population')
        plt.legend()
        plt.grid(True)

        # Travelers from B to A
        plt.subplot(4, 1, 4)
        plt.plot(range(weeks + 1), Sba, 'b-', label='Susceptible in A from B')
        plt.plot(range(weeks + 1), Iba, 'r-', label='Infected in A from B')
        plt.plot(range(weeks + 1), Rba, 'g-', label='Recovered in A from B')
        plt.title('Travelers from B to A')
        plt.xlabel('Weeks')
        plt.ylabel('Population')
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.savefig('SIR_Model_Travelers.png')
        plt.show()


if __name__ == "__main__":
    # Example parameters
    Na = 1000  # Population of group A
    Nb = 1000  # Population of group B
    Ia0 = 5  # Initial infected in group A
    Ib0 = 0  # Initial infected in group B
    weeks = 50  # Simulation duration in weeks

    Reca = 0.1  # Recovery rate for group A
    Recb = 0.1  # Recovery rate for group B
    Ka = 0.3  # Transmission coefficient for group A
    Kb = 0.3  # Transmission coefficient for group B
    Mab = 0.1  # Proportion of A traveling to B
    Mba = 0.1  # Proportion of B traveling to A

    SIR_Model(Na, Nb, Ia0, Ib0, Reca, Recb, Ka, Kb, Mab, Mba, weeks)