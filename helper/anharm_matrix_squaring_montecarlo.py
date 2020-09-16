import math, numpy, pylab
from numpy import *

# Define the anharmonic (quartic) potential
def V_anharmonic(x, gamma, kappa):
    V = x**2 / 2 + gamma * x**3 + kappa * x**4
    return V

# Free off-diagonal density matrix:
def rho_free(x, xp, beta):
    return (math.exp(-(x - xp) ** 2 / (2.0 * beta)) /
            math.sqrt(2.0 * math.pi * beta))

# Harmonic density matrix in the Trotter approximation (returns the full matrix):
def rho_anharmonic_trotter(grid, beta):
    return numpy.array([[rho_free(x, xp, beta) * \
            numpy.exp(-0.5 * beta * (V_anharmonic(x, -g, g) + V_anharmonic(xp, -g, g))) \
            for x in grid] for xp in grid])

# Exact harmonic oscillator quantum position distribution:
def p_quant(x, beta):
    p_q = sqrt(tanh(beta / 2.0) / pi) * exp(- x**2.0 * tanh(beta / 2.0))
    return p_q

# Perturbative energy levels
def Energy_pert(n, cubic, quartic):
    return n + 0.5 - 15.0 / 4.0 * cubic **2 * (n ** 2 + n + 11.0 / 30.0) \
         + 3.0 / 2.0 * quartic * (n ** 2 + n + 1.0 / 2.0)

# Partition function obtained using perturbative energies
def Z_pert(cubic, quartic, beta, n_max):
    Z = sum(math.exp(-beta * Energy_pert(n, cubic, quartic)) for n in range(n_max + 1))
    return Z

# Construct the position grid:
x_max = 5 #maximum position value
nx = 100 #number of elements on the x grid
dx = 2.0 * x_max / (nx - 1) #position differential
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)] #position grid

beta_tmp = 2.0 ** (-5) # initial (low) value of beta (power of 2) (high temperature)
beta = 2.0 ** 1 # actual value of beta (power of 2)

#g = 1.0 #-cubic and quartic coefficient

for g in [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]:
    
    Z_p = Z_pert(-g, g, beta, 15)
    
    rho = rho_anharmonic_trotter(x, beta_tmp) # density matrix at initial (low) beta (Trotter decomp.)

    # Reduce the temperature by the convolution property (matrix squaring):
    while beta_tmp < beta:
        rho = numpy.dot(rho, rho) #matrix squaring (convolution)
        rho *= dx #also multiply by the differential since we are in position representation
        beta_tmp *= 2.0 #reduce the temperature by a factor of 2
        #print 'beta: %s -> %s' % (beta_tmp / 2.0, beta_tmp)

    # Output position distribution pi(x) at the final beta onto a file:
    Z = sum(rho[j, j] for j in range(nx + 1)) * dx #partition function
    pi_of_x = [rho[j, j] / Z for j in range(nx + 1)] #the diagonal element of the density matrix
    f = open('data_anharm_matrixsquaring_beta' + str(beta) + '.dat', 'w')
    for j in range(nx + 1):
        f.write(str(x[j]) + ' ' + str(rho[j, j] / Z) + '\n')
    f.close()

    # Plot the obtained final position distribution:
    T = 1 / beta
    x = linspace(-x_max, x_max, nx+1)
    y2 = [V_anharmonic(a, -g, g) for a in x]
    y1 = [p_quant(a, beta) for a in x]
    pylab.plot(x, y2, c='gray', linewidth=2.0, label='Anharmonic potential')
    pylab.plot(x, y1, c='green', linewidth=2.0, label='Harmonic exact quantum')
    pylab.plot(x, pi_of_x, c='red', linewidth=4.0, label='Anharmonic matrix squaring')
    pylab.ylim(0,1)
    pylab.xlim(-2,2)
    pylab.title('Anharmonic oscillator position distribution at $T=$%.2f' % T, fontsize = 13)
    pylab.xlabel('$x$', fontsize = 15)
    pylab.ylabel('$\pi(x)$', fontsize = 15)
    pylab.legend()
    pylab.grid()
    pylab.savefig('plot_T_%.2f_anharm_g_%.1f_prob_matrix_squaring.png' % (T,g))
    pylab.show()
    print 'g =', g, 'Perturbative partition function:', Z_p, 'Monte Carlo partition function', Z