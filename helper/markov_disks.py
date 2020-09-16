import random, math, cmath, pylab, os

#Introduce the periodic boundary conditions via the modular distance between two 2D vectors x, y:
def dist(x,y):
    d_x = abs(x[0] - y[0]) % 1.0 #distance between the first compononents of two vectors in modulo 1
    d_x = min(d_x, 1.0 - d_x) #the modular distance is the minimum of cases due to periodicity
    d_y = abs(x[1] - y[1]) % 1.0 #distance between the second compononents of two vectors in modulo 1
    d_y = min(d_y, 1.0 - d_y)
    return  math.sqrt(d_x**2 + d_y**2) #returns the modular Euclidian distance magnitude
def delx_dely(x, y):
    d_x = (x[0] - y[0]) % 1.0
    if d_x > 0.5: d_x -= 1.0
    d_y = (x[1] - y[1]) % 1.0
    if d_y > 0.5: d_y -= 1.0
    return d_x, d_y #returns the modular distance vector

#Define the integer square root function (to use in range() function)
def isqrt(n):
    i = int(math.sqrt(n) + 0.5)
    if i**2 == n:
        return i
    raise ValueError('input was not a perfect square')

def Psi_6(L, eta):
    sigma = math.sqrt(eta / N / math.pi) #radius
    sum_vector = 0j #initialise
    for i in range(N):
        vector  = 0j #initialise
        n_neighbor = 0 #initialise the number of neighbours
        for j in range(N):
            if dist(L[i], L[j]) < 2.8 * sigma and i != j: #neighbour condition
                n_neighbor += 1 #increase the number of neighbours if condition satisfied
                dx, dy = delx_dely(L[j], L[i]) #find the distance between neighbouring disks
#find the angle between x axis and vector (dx,dy) by treating it as a complex number and finding its phase
                angle = cmath.phase(complex(dx, dy)) 
                vector += cmath.exp(6.0j * angle) #the local orientational order parameter
        if n_neighbor > 0:
            vector /= n_neighbor 
        sum_vector += vector #the global order parameter (sum of local order parameters)
    return sum_vector / float(N)

#markov-chain sampling for 4 hard disks with periodic boundary conditions
def markov_disks(L, delta, eta):
    sigma = math.sqrt(eta / N / math.pi) #radius
    condition = True #initialise the condition of overlap as True
    while condition == True:
        a = random.choice(L)
        b = [(a[0] + random.uniform(-delta, delta)) % 1.0, (a[1] + random.uniform(-delta, delta)) % 1.0]
        min_dist = min(dist(b, c) for c in L if c != a)
        if not min_dist < 2.0 * sigma:
            a[:] = b
            condition = False #overlap condition is false, i.e. there is no overlap
            break
    return L

#Inputs of the markov_disks_box function:
N = 64 #number of disks
eta = 0.72 #disk density
delta = 0.5 * math.sqrt(eta / N / math.pi)
#Obtain the initial disk configuration L from a file if it exists; otherwise generate it.
filename = 'disk_configuration_N%i_eta%.2f.txt' % (N, eta)
if os.path.isfile(filename): #Returns True if "filename" is an existing regular file.
    f = open(filename, 'r')
    L = []
    for line in f:
        a, b = line.split()
        L.append([float(a), float(b)])
    f.close()
    print 'starting from file', filename
else:
    delxy = 1 / math.sqrt(N) / 2
    L = [[delxy + i * 2 * delxy, delxy + j * 2 * delxy] for i in range(isqrt(N)) for j in range(isqrt(N))]
    print 'starting from a new configuration'
    
n_steps = 10000 #number of steps

#Figure output:
def show_conf(L, eta, title, fname):
    sigma = math.sqrt(eta / N / math.pi) #radius
    pylab.axes()
    for [x, y] in L:
        #the following loops generate 6 disks on a square grid to implement periodic boundary conditions
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                cir = pylab.Circle((x + ix, y + iy), radius=sigma,  fc='r')
                pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    pylab.title(title)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    #pylab.savefig(fname)
    pylab.show()
    pylab.close()   
    
del_eta = 0.02 #the amount the density is reduced after each run

condition_eta = False #the condition that eta is greater than or equal to 0.2
run = 0  
Psi = []

while condition_eta == False:
    n_avg = 0
    Psi_avg = 0
    run += 1
    #draw the intermediate configuration
    #show_conf(L, eta, 'Intermediate configuration N=%i eta=%.2f' % (N, eta), 'filename.png') 
    for step in range(n_steps):
        #pylab.clf()
        #show_conf(L, sigma, 'Sample Configuration', 'disks.png')
        markov_disks(L, delta, eta)    
        #Calculate the absolute value of the complex order parameter at each 100 steps:
        if round(step / 100) == step / 100:
            n_avg += 1
            Psi_avg += abs(Psi_6(L, eta))
    Psi.append(Psi_avg / n_avg) 
    print Psi_avg / n_avg
    
    eta -= del_eta #Reduce the density by a small amount after each run;
    print run, eta #print the current run number
    if eta <= 0.2: #if eta is less than 0.2, then break
        conditon_eta = True
        break
        
    #Write the final configuration on a file:
    filename = 'disk_configuration_N%i_eta%.2f.txt' % (N, eta)
    f = open(filename, 'w')
    for a in L:
       f.write(str(a[0]) + ' ' + str(a[1]) + '\n')
    f.close()

#Plot of the average of |Psi_6| as function of the density eta
x = pylab.linspace(0.72, 0.2, len(Psi))
pylab.plot(x, Psi, 'black', linewidth=2.0)
pylab.xlabel('$\eta$', fontsize=14)
pylab.ylabel('$\Psi_6(\eta)$', fontsize=14)
pylab.title('The global order parameter for hard disks')
pylab.savefig('global-order_parameter_disks.png')
pylab.show()