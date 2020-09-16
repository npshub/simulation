import math, random, pylab, os

def show_path(x, step, title):
    fig = pylab.figure(figsize=(6, 10))
    pylab.clf()
    path = x[:] #Final position is the same as the initial position. 
    #Note that this notation appends the first element of x as a new element to x 
    y_axis = range(len(x)) #construct the imaginary time axis
    #Plot the paths
    pylab.plot(path, y_axis, 'ko-')
    pylab.xlabel('$x$', fontsize = 14)
    pylab.ylabel('$\\tau$', fontsize=14)
    pylab.title(title, fontsize=14)
    pylab.grid()
    pylab.savefig('levy_snapshot_%i.png' % step)
    
def V_harm(x):
    pot = x**2 / 2
    return pot

def V_anharm(x, cubic, quartic):
    pot = x**2 / 2 + cubic * x**3 + quartic * x**4
    return pot

def levy_harmonic_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        Ups1 = 1.0 / math.tanh(dtau) + \
               1.0 / math.tanh(dtau_prime)
        Ups2 = x[k - 1] / math.sinh(dtau) + \
               xend / math.sinh(dtau_prime)
        x.append(random.gauss(Ups2 / Ups1, \
               1.0 / math.sqrt(Ups1)))
    return x

beta = 20.0
N = 100
dtau = beta / N
delta = 1.0
n_steps = 10000
x = [1.0] * N
cubic, quartic = -1, 1
data = []
for step in range(n_steps):
    x_old = x[:]
    Ncut = 11 #pick a prime number ?
    x = x[Ncut:] + x[:Ncut] #slices the x data from Ncut'th element and swaps the two parts, hence x[0] changes
    x_new = levy_harmonic_path(x[0], x[Ncut], dtau, Ncut) + x[Ncut:] #sample x points using harmonic Levy sampling
    #This is why the free and harmonic density matrix part is not necessary in the Metropolis acceptance rate
    #Metropolis algorithm
    if random.uniform(0,1) < math.exp(sum(-(V_anharm(a, cubic, quartic) - V_harm(a)) * dtau for a in x_new)) \
                            / math.exp(sum(-(V_anharm(a, cubic, quartic) - V_harm(a)) * dtau for a in x_old)): 
        x = x_new[:]
    k = random.randint(0, N - 1) 
    data.append(x[k])
        
#Write the final path configuration on a file:
filename = 'improved_levy_anharmonic_path_configuration_levy_N%i.txt' % N
f = open(filename, 'w')
for i in range(N):
    f.write(str(x[i])+ '\n')
f.close()

show_path(x, step, 'Levy anharmonic final configuration, N = %i' %N)
pylab.show()

pylab.clf()
fig = pylab.figure(figsize=(7, 5))
T = 1 / beta
list_x = [0.1 * a for a in range (-30, 31)]
y = [V(a, cubic, quartic) for a in list_x]
pylab.plot(list_x, y, c='gray', linewidth=4.0, label='Anharmonic potential')
pylab.hist(data, normed=True, bins=100, label='QMC')
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('levy_anharmonic_path_improved (beta=%s, N=%i, cubic=%.1f, quartic=%.1f)' % (beta, N, cubic, quartic))
pylab.xlim(-2, 2)
pylab.ylim(0, 0.9)
pylab.savefig('plot_C2_beta%s.png' % beta)
pylab.show()