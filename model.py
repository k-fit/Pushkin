import math
import sys
from pyOpt import Optimization
from pyOpt import SLSQP
from pyOpt import SOLVOPT 
from scipy import special
from collections import Counter


# This is an implementation of the BG-NBD model as described in the paper "Counting Your Customers the Easy Way" by Fader, Hardie, and Lee, 2005.
# The data is a sample set, provided by them, supposedly corresponding to CD-NOW data.
# There are only three pieces of data per customer.  x = number of transactions observed in the time period [0, T] and tx, the time since the last transaction.  
# The maximum value for T is 39 weeks.  
# Time is, confusingly, measured in weeks.

FILENAME = 'bg_data.dat'

# This is the log likelihood of observing the data given the model parameters
def objfun(x):
    r = x[0]
    alpha = x[1]
    a = x[2]
    b = x[3]
    
    
    data = open(FILENAME, 'r')

    LL = 0.0
    g = []
    g = [0.0]*4
    g[0] = -1*(x[0] - .00001)
    g[1] = -1*(x[1] - .00001)
    g[2] = -1*(x[2] - .00001)
    g[3] = -1*(x[3] - .00001)
    #g[0] *= -1
    #g[1] *= -1
    #g[2] *= -1
    #g[3] *= -1

    for line in data:
        values = line.split('\t')
        if values[0] == 'ID': continue
    
        values = [v.strip() for v in values]
    
        x = float(values[1])
        tx = float(values[2])
        T = float(values[3])
    
        #if alpha > 0:
        log_A1 = math.lgamma(r + x) - math.lgamma(r) + r * math.log(alpha)
        #else:
        #log_A1 = -10000
        log_A2 = math.lgamma(a + b) + math.lgamma(b + x) - math.lgamma(b) - math.lgamma(a + b + x)

        log_A3 = -1 * (r + x) * math.log(alpha + T)

        if x > 0: 
            log_A4 = math.log(a) - math.log(b + x - 1) - (r + x) * math.log(alpha + tx)
            d = 1
        else:
            log_A4 = 0
            d = 0
    
        try:
            contribution = log_A1 + log_A2 + math.log(math.exp(log_A3) + d * math.exp(log_A4))
        except:
            print log_A1, log_A2, math.exp(log_A3), d, math.exp(log_A4)
            sys.exit(1)

        LL += contribution
    fail = 0
    
    # pyOpt attempts to minimize a function, so we multiply the result by -1
    return  -1*LL, g, fail

def runoptimizer():
    opt_prob = Optimization('TP37 Constrained Problem',objfun)
    opt_prob.addObj('LL')
    opt_prob.addVar('x1','c',lower=0.01,upper=10.0,value=1.0)
    opt_prob.addVar('x2','c',lower=0.01,upper=10.0,value=1.0)
    opt_prob.addVar('x3','c',lower=0.01,upper=10.0,value=1.0)
    opt_prob.addVar('x4','c',lower=0.01,upper=10.0,value=1.0)

    opt_prob.addConGroup('g', 4, 'i')

    # sanity check
    print opt_prob
    print objfun([1.0,1.0,1.0,1.0])

    # other optimization methods can be used here - we use sequential least squares programming
    slsqp = SLSQP() 
    [fstr, xstr, inform] = slsqp(opt_prob)

    print opt_prob.solution(0)
    return [v.value for v in opt_prob.solution(0).getVarSet().values()]


# Given the model parameters we can compute the expected number of transactions for a random customer in a time period of length tt
# These expected values are stored in a dict, indexed by time
def getExpectations(r, alpha, a, b, numdays):
    expectations = dict()
    for tt in range(1, numdays + 1):
        t = tt / 7.0
        x = special.hyp2f1(r, b, a + b - 1, t / (alpha + t))
        expectations[tt] = ((a + b - 1) / (a - 1)) * ( 1 - ( (alpha / (alpha + t)) ** r ) * x)
    return expectations

# This is a frequency table of the number of users who made their first purchase on day n.
def getTimeDistributions():
    file = open(FILENAME, 'r')
    X = list()

    for line in file:
        values = line.split('\t')
        if values[0] == 'ID': continue
    
        values = [v.strip() for v in values]
    
        x = float(values[1])
        tx = float(values[2])
        T = float(values[3])
       
        numdays = round((39 - T) * 7)

        X.append(int(numdays))

    return Counter(X)

# Rather than compute the expected transactions for a randomly chosen individual,
# We would rather estimate and predict the number of transactions made by the entire popultation we have in our data set.
# The answer is again stored in a dict object, indexed by days.
def getTransactions(firstPurchaseCounts, expectations, numdays):
    total = dict()
    for tt in range(1, numdays + 1):
        t = tt / 7.0
        s = 1
        tot = 0
        while (tt > s):
            tot += firstPurchaseCounts[s ] * expectations[tt - s]
            s += 1
        total[t] = tot
    return total

# This calculation is separate from those above.  For a particular customer, we can predict his/her expected purchases using the past transaction information we have.
# Given information of purchases in time [0, T] (ie x, and tx) and the model parameters, we can compute expected purchases in the next t weeks.
def conditionalExpectation(t, x, tx, T, r, alpha, a, b):
    E = (a + b + x - 1) / (a - 1)

    E *= 1 - (( (alpha + T) / (alpha + T + t) ) ** (r + x)) * special.hyp2f1(r + x, b + x, a + b + x - 1, t / (alpha + T + t))

    if x > 0:
        E /= 1 + (a / (b + x - 1)) * (((alpha + T) / (alpha + tx))**(r + x))

    return E

def main():
    # Solve for model parameters
    [r, alpha, a, b] = runoptimizer()
    #r, alpha, a, b = 0.24259608786, 4.41356726546, 0.792913903499, 2.42599309994;

    print "Model Parameters (r, alpha, a, b): ",
    print r, alpha, a, b
    
    # predict out to week 78 - nothing special about 78
    numdays = 7 * 78

    # get expected value for a randomly chosen customer, indexed by time
    expectations = getExpectations(r, alpha, a, b, numdays)
    
    # log the frequency of customers who started on a given day
    firstPurchaseCounts = getTimeDistributions()

    # use the above two computations to compute the expected purchases of our population
    expectedTransactions = getTransactions(firstPurchaseCounts, expectations, numdays)

    # TODO : visualize expected transactions with prediction against actual transactions
    
    # We can also compute expected value per customer.  Here it's just done for the first customer.
    f = open('data.txt', 'r')
    f.readline()
    values = f.readline().strip().split('\t')
    x = int(values[1])
    tx = float(values[2])
    T = float(values[3])
    t = round(78 - T)
    print "Expected future purchases for customer 1: ", 
    print conditionalExpectation(t, x, tx, T, r, alpha, a, b)


if __name__ == '__main__':
    main()
