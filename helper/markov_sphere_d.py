import pylab, random, math, sys, os

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

#calculate the analytical values:
def V_sph(dim):
    return math.pi ** (dim / 2.0) / math.gamma(dim / 2.0 + 1.0)
def Q(dim):
    return V_sph(dim) / V_sph(dim-1)

#monte carlo algorithm:
def markov_sphere(d, n_runs, n_trials):
    x = [0.0] * d #initial point
    delta = 0.1
    Q_avg_d = 0 #initialise the average Q(d)
    old_radius_square = sum(i**2 for i in x)
    data = []
    for j in range(n_runs):
        n_hits = 0
        if d == 1:
            Q_avg_d = 1 #multiplicative identity (not the actual value)
            break
        for i in range(n_trials):
            k = random.randint(0, d - 2)
            x_old_k = x[k]
            x_new_k = x_old_k + random.uniform(-delta, delta)  
            new_radius_square = old_radius_square + x_new_k ** 2 - x_old_k ** 2
            if new_radius_square < 1.0: #check whether the position is inside the unit (d-1)-sphere
                old_radius_square = new_radius_square #update the radius
                x[k] = x_new_k #update a component
                x[d-1] = random.uniform(-1, 1) #sample the d'th coordinate (inside the unit cylinder)                
            if old_radius_square + x[d-1]**2 < 1.0: #is x inside the unit d-sphere?
                n_hits += 1 #raise the number of hits in case of accepted sample                                
                data.append(math.sqrt(old_radius_square  + x[d-1]**2)) #generate the radius histogram data
                
        Q_avg_d += 2.0 * n_hits / float(n_trials) / n_runs #take the average of Q's
        
    print 'Analytical: Q(d=%i)' % d, Q(d)
    if d != 1:
        print 'Numerical: <Q(d=%i)> =' % d, Q_avg_d 
    else:
        print 'There is no cylinder in 1 dimension!'
    return Q_avg_d, data

#Calculate the volume of the d_max dimensional sphere through the formula V(n)=Q(n)...Q(2)V(2):
def V_sph_markov(dim, n_trials, n_runs):
    V = []
    V_d = V_sph(1) #initialise the volume of "d_max" dimensional sphere
    print 'The analytical value for the volume of the unit %i-sphere is =' % dim, V_sph(dim), '.'
    print '___________________________________________'
    for d in range(1, dim + 1):    
        Q_avg_d, data = markov_sphere(d, n_runs, n_trials)
        V_d *= Q_avg_d
        V.append(V_d) #save the volume for each d value along the way until d_max
        print 'Analytical: V_sph(d=%i)' % d, V_sph(d)
        print 'Numerical: V_sph(d=%i)' % d, V_d
        print '___________________________________________'
    print 'After %i runs' % n_runs, 'each consisting of %i trials' % n_trials, ', the numerical result for the volume of the unit %i-sphere (with analytical value' % dim, V_sph(dim), ') is found to be', V_d, '.'
    return V

d_max = 200 #maximum dimension
n_trials = 1000 #each dimension iteration takes "n_trials" number of iterations
n_runs = 500
V = V_sph_markov(d_max, n_trials, n_runs)

#Error calculation:
d = 20 #dimension of the sphere
n_runs = 10
Error = []
trials = [10, 100, 1000, 10000, 100000]

for trial in trials:
    V_avg, V_avg_square = 0.0, 0.0
    for run in range(n_runs):
        with HiddenPrints():
            v = V_sph_markov(d, 1, trial)
        V_avg += v[d-1] / n_runs
        V_avg_square += v[d-1] ** 2.0 / n_runs
    print V_avg
    Error.append(math.sqrt(V_avg_square - V_avg ** 2) / pow(n_runs, 0.5))   
print Error
Dif = [] * len(trials)
for i in range(len(trials)):
    Dif.append(abs(Error[i] - V_sph(d)))
print Dif