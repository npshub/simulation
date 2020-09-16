import random

def markov_disks_box(L, delta, sigma):
    condition = True
    while condition == True:
        a = random.choice(L)
        b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
        min_dist = min((b[0] - c[0]) ** 2 + (b[1] - c[1]) ** 2 for c in L if c != a)
        box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma
        if not (box_cond or min_dist < 4.0 * sigma ** 2):
            a[:] = b
            condition = False
            break
    return L

#inputs of the markov_disks_box function:
#initial positions of the disks to startup the Markov-chain
L = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]] 
delta = 0.1
sigma = 0.15 #radius

n_steps = 1000000
del_xy = 0.05 #"uncertainty"

#Define the "marked" configurations:
conf_a = ((0.30, 0.30), (0.30, 0.70), (0.70, 0.30), (0.70,0.70))
conf_b = ((0.20, 0.20), (0.20, 0.80), (0.75, 0.25), (0.75,0.75))
conf_c = ((0.30, 0.20), (0.30, 0.80), (0.70, 0.20), (0.70,0.70))
configurations = [conf_a, conf_b, conf_c] #list the configurations

hits = {conf_a: 0, conf_b: 0, conf_c: 0} #initialise the number of times each marked configuration occurs

for run in range(n_steps):
    x_vec = direct_disks_box(4, sigma) #generates a random sample by direct sampling
    for conf in configurations: #run a loop iterating over the given 3 configurations
    #condition that a randomly generated configuration L is the same as a, b or c up to uncertainty of del_xy
        condition_hit = True
        for b in conf: #run a loop iterating over the 4 disk coordinates in a specific configuration
            #If the max(x distance and y distance between a disk in L and a disk in conf_a,b,c) 
            #is less than the given uncertainty del_xy, then we treat the two disks as in the same location.
            #Note that the "any two disks" condition is realised by minimising over all 4 disks in a 
            #randomly sampled configuration L.
            condition_b = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a in x_vec) < del_xy  
            #The following logical variable is 1 only if there exists 4 disk pairs are within del_xy range.
            #If at least any one of the disks does not have a pair within del_xy, then it is 0.
            condition_hit *= condition_b #multiplies condition_b's (for all 4 disks)
        #If the current L and a, b or c are the same up to uncertainty del_xy, then increase:
        if condition_hit:
            hits[conf] += 1

for conf in configurations:
    print conf, hits[conf] #Print the configurations and the number of times they occured.