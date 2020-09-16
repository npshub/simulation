%pylab inline
import math, random, pylab

# Define the anharmonic (quartic) potential
def V_anharmonic(x, gamma, kappa):
    V = x**2 / 2 + gamma * x**3 + kappa * x**4
    return V

def rho_free(x, y, beta): # free off-diagonal density matrix
    return math.exp(-(x - y) ** 2 / (2.0 * beta)) 

def read_file(filename):
    list_x = []
    list_y = []
    with open(filename) as f:
        for line in f:
            x, y = line.split()
            list_x.append(float(x))
            list_y.append(float(y))
    f.close()
    return list_x, list_y

beta = 4.0
g = 1.0 #-cubic and quartic coefficients

T = 1 / beta
N = 16                                            # number of imaginary times slices
dtau = beta / N
delta = 1.0                                       # maximum displacement on one slice
n_steps = 1000000                                 # number of Monte Carlo steps
x = [0.0] * N                                     # initial path
hist_data = []
for step in range(n_steps):
    k = random.randint(0, N - 1)                  # random slice
    knext, kprev = (k + 1) % N, (k - 1) % N       # next/previous slices
    x_new = x[k] + random.uniform(-delta, delta)  # new position at slice k
    old_weight  = (rho_free(x[knext], x[k], dtau) *
                   rho_free(x[k], x[kprev], dtau) *
                   math.exp(-dtau * V_anharmonic(x[k], -g, g)))
    new_weight  = (rho_free(x[knext], x_new, dtau) *
                   rho_free(x_new, x[kprev], dtau) *
                   math.exp(-dtau * V_anharmonic(x_new ,-g, g)))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x[k] = x_new
    if step % 10 == 0:
        hist_data.append(x[0])
        
# Figure output:
list_x, list_y = read_file('data_anharm_matrixsquaring_beta' + str(beta) + '.dat')
v = [V_anharmonic(a, -g, g) for a in list_x]
pylab.plot(list_x, v, c='gray', linewidth=2.0, label='Anharmonic potential')
pylab.plot(list_x, list_y, c='red', linewidth=4.0, label='path integral Monte Carlo')
pylab.hist(hist_data, 100, normed = 'True', label='matrix squaring') #histogram of the sample
pylab.ylim(0,1)
pylab.xlim(-2,2)
pylab.title('Position distribution at $T=%.2f$, $\gamma_{cubic}=%.2f$, $\gamma_{quartic}=%.2f$' % (T,-g,g), fontsize = 13)
pylab.xlim(-2.0, 2.0) #restrict the range over which the histogram is shown
pylab.xlabel('$x$', fontsize = 15)
pylab.ylabel('$\pi(x)$', fontsize = 15)
pylab.legend()
pylab.savefig('plot_T_%.2f_anharm_g_%.1f_prob_path_int.png' % (T,g))
pylab.show()