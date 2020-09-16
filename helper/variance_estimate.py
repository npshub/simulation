import random, math
n_trials = 400000
n_hits = 0
var = 0.0
var_est = 0.0
Obs_exp = 0.0
Obs2_exp = 0.0
for iter in range(n_trials):
    x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
    Obs = 0.0
    if x**2 + y**2 < 1.0:
        n_hits += 1
        Obs = 4.0
    Obs_exp+=Obs/n_trials #expectation value of the observable
    Obs2_exp+=(Obs)**2/n_trials #expectation value of the square of the observable
    var+=(Obs-math.pi)**2/n_trials #calculation of the variance if the value of pi had been known beforehand    
var_est = Obs2_exp - Obs_exp**2 #best estimate to the variance
print 4.0 * n_hits / float(n_trials), math.sqrt(var), math.sqrt(var_est)