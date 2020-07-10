import ecos
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt, sin
from scipy.sparse import csc_matrix
from numpy import linalg, concatenate

class PDG:
    def __init__(self, m_f, rho_1, rho_2, alpha, gamma, theta):
        # vessel specific constant values
        self.m_f   = m_f
        self.rho_1 = rho_1
        self.rho_2 = rho_2
        self.alpha = alpha
        self.gamma = gamma
        self.theta = theta
        self.B     = concatenate((np.zeros((3, 3)), np.identity(3), np.zeros((1, 3))), axis=0)
        self.B     = concatenate((self.B, np.array([[0, 0, 0, 0, 0, 0, -self.alpha]]).T), axis=1)
        self.zeta  = 1 - .5 * (3 - sqrt(5))

        # user-defined hyperparameters
        self.waitTime = 0
        self.predTime = 0

    def INIT_CONSTANTS(self, g, w):
        # takes planet specifc values, gravity and planet rotation rate, as inputs and initializes
        # matrices used in equations of motions that use these planet-specific values
        self.g = g

        s_W = np.array([
                        [w[2]**2+w[1]**2,   -w[1]*w[0],         -w[2]* w[0]     ,   0,          2*w[2], -2*w[1]],
                        [-w[1]*w[0],        w[2]**2+w[0]**2,    -w[2]*w[1]      ,   -2*w[2],    0,       2*w[0]],
                        [-w[2]*w[0],        -w[2]*w[1],         w[1]**2+w[0]**2 ,   2*w[1],     -2*w[0],      0],
                       ])

        self.A = concatenate((np.zeros((3, 3)), np.identity(3)), axis=1)
        self.A = concatenate((self.A, s_W, np.zeros((1, 6))))
        self.A = concatenate((self.A, np.zeros((7, 1))), axis=1)

    def PREDICT_STATE_VECTOR(self, predTime, position, velocity, mass, currentEta):
        # updates current state vector
        if predTime != self.predTime:
            self.phi = self.mat_exp(self.dt)
            self.psi = np.trapz([np.dot(self.mat_exp(.02 * tau), self.B) for tau in range(50)],
                                axis=0, dx=.02 * predTime)
            self.predTime = predTime

        x_0   = np.array([position[1], position[0], position[2], 
                          velocity[1], velocity[0], velocity[2], 
                          log(mass)]) 

        self.x  = np.dot(self.phi, x_0) + np.dot(self.psi, currentEta + np.array([self.g[7], 0, 0, 0]))
        self.x  = np.concatenate((self.x, np.array(currentEta)))
        self.m0 = exp(self.x[6])

    def SCHEDULE_PDG(self):
        '''
         Optimizes time to ignition, solvetime, and final landing error bound
            Parameters:
                x0 (array):         initial state vector at time of call
            Returns:
                waitTime (double):  optimal time to ignition, in seconds
                tSolve (double):    optimal solvetime, in seconds
                dMax (double):      optimal maximum landing error
        '''

        tLow, tHigh = 0, sqrt(2 * self.x[0] / -g[7]) / self.dt

        waitTime_1 = int(tHigh - self.zeta * (tHigh - tLow))
        waitTime_2 = int(tLow + self.zeta * (tHigh - tLow))
        fFuel_1 = self.MIN_DISTANCE(waitTime_1)
        fFuel_2 = self.MIN_DISTANCE(waitTime_2)
        
        while abs(waitTime_1 - waitTime_2) > int(1 / self.dt):       
            if fFuel_1 == None or (fFuel_2 != None and fFuel_1 > fFuel_2):
                tLow = waitTime_1
                waitTime_1, fFuel_1 = waitTime_2, fFuel_2
                
                waitTime_2 = int(tLow + self.zeta * (tHigh - tLow))
                fFuel_2 = self.MIN_DISTANCE(waitTime_2, schedule_pdg=False)

            else:
                tHigh = waitTime_2
                waitTime_2, fFuel_2 = waitTime_1, fFuel_1
                
                waitTime_1 = int(tHigh - self.zeta * (tHigh - tLow))
                fFuel_1 = self.MIN_DISTANCE(waitTime_1, schedule_pdg=False)
        
        self.waitTime = waitTime_2 * self.dt  

    def MIN_DISTANCE(self, waitTime=0, schedule_pdg=False):

        ''' Optimizes solvetime and final landing distance.
            Parameters:
                waitTime (double): time to ignition, in seconds
            Returns:
                t1 (double):    optimal solvetime, in seconds
                dOpt (double):  optimal final landing error
        '''

        # creates equations of motion used in equality matrix, if waitTime is greater
        # than 0, finds state vector after waitTime time passes before ignition
        
        self.create_equations_of_motion()

        x_s, eta = self.x, self.x[7:11]
        for _ in range(int(waitTime)):
            x_s = concatenate(( np.dot(self.omega, (x_s + self.g * self.dt)), eta ))
        m_s = np.exp(x_s[6])                                  

        # upper bound for line search is time it takes to run out of fuel, lower bound is zero

        tHigh = ((m_s - self.m_f) / (self.alpha * self.rho_2 * self.dt)) 
        tLow  = 0

        t1 = int(tHigh - self.zeta * (tHigh - tLow))
        t2 = int(tLow + self.zeta * (tHigh - tLow))
        err_1, _, fDist_1 = self.run_ecos(t1, x_s, m_s)
        err_2, _, fDist_2 = self.run_ecos(t2, x_s, m_s)

        # conducts golden search until t1 and t2 both yield optimal solutions, and the final 
        # landing distance that corresponds to both t1 and t2 are close enough to each other, 
        # returns t1
    
        while err_1 != 0 or err_2 != 0 or abs(fDist_1 - fDist_2) > .01:
            # if t1 throws an error or t2 is valid and distance 1 is greater than distance two, 
            # move the lower time bound to the right, otherwise, move upper time bound to the left
            if err_1 != 0 or (err_2 == 0 and fDist_1 > fDist_2):
                tLow, t1          = t1, t2
                err_1, fDist_1    = err_2, fDist_2 
                t2                = int(tLow + self.zeta * (tHigh - tLow))               
                err_2, _, fDist_2 = self.run_ecos(t2, x_s, m_s) 

            else:
                tHigh, t2         = t2, t1
                err_2, fDist_2    = err_1, fDist_1
                t1                = int(tHigh - self.zeta * (tHigh - tLow))
                err_1, _, fDist_1 = self.run_ecos(t1, x_s, m_s)

            # if the line search found no feasible solutions, return None
            if err_1 != 0 and err_2 != 0 and abs(t1 - t2) < int(1/self.dt):
                print(t1, t2)
                return None

        if schedule_pdg: return fDist_1
        else:  
            #print("Maximum distance = ", fDist_1)
            self.dMax = fDist_1          
            return t1 * self.dt

    def MIN_FUEL(self, tSolve):
        # creates equality and inequality matrices used in ecos, calls ecos with negative final
        # fuel is the minimization function, and returns a list of eta vectors to main_control

        self.create_equations_of_motion()

        G, h, q, l, e = self.soc_constraints(int(tSolve/self.dt), self.m0)
        A_mat, b      = self.equality_constraints(int(tSolve/self.dt), self.x)
        c             = np.zeros(11 * int(tSolve/self.dt) + 1)

        c[11*(int(tSolve/self.dt)-1) + 6] = -1
        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                                verbose=False, abstol=1e-4, feastol=1e-4, reltol=1e-4)
        sol = solution['x']      
        
        # from the optimized solution found by ecos, a list is 
        # created of the thrust vectors at every time step.
        
        eta = []                                                    
        
        for t in range(int(tSolve/self.dt)):                                   
            thrustVect = sol[11*t + 7:11*t + 10]                       
            thrustMagn = linalg.norm(thrustVect)
            
            eta.extend(thrustVect)
            eta.append(thrustMagn)
            
        return sol, eta

    def create_equations_of_motion(self):
        #  Define matrices used in ODE system for every timestep; used to construct equality matrix.

        phi   = self.mat_exp(self.dt)
        psi   = np.trapz([np.dot(self.mat_exp(.002 * tau * self.dt), self.B) for tau in range(500)],
                            axis=0, dx=.002 * self.dt)
        omega = concatenate((phi, psi), axis=1)
        E     = concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)

        # rewrites matrices in csc format
        
        self.rowList, self.colList, self.valList = [], [], []
        for r in range(len(E)):
            for c in range(len(E[0])):
                if E[r, c] != 0:
                    self.rowList.append(r)
                    self.colList.append(c)
                    self.valList.append(E[r, c])
        
        self.bVect = np.dot(omega, self.g)

    def mat_exp(self, x):

        '''  Taylor polynomial approximation of matrix exponential with scalar multiplication
        Parameters:
            A (matrix): array
            x (double): scalar
        Returns:
            expMat:     approximate matrix exponential as numpy array 
        '''
        
        expMat = np.zeros_like(self.A)
        for i in range(30):
            expMat = expMat + np.linalg.matrix_power(self.A, i) * x**i / factorial(i)
        return expMat

    def run_ecos(self, tSolve, x_s, m_s):
        minDistance   = True
        G, h, q, l, e = self.soc_constraints(int(tSolve/self.dt), m_s, minDistance=True)
        A_mat, b      = self.equality_constraints(int(tSolve/self.dt), x_s)
        c             = np.zeros(11 * int(tSolve/self.dt) + 1)
        
        #  Optimize slack variable; final distance is less than or equal to slack variable. 
        #c[-1] = 1
        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                              verbose=False, abstol=1e-4, feastol=1e-4, reltol=1e-4)
        fDist = np.linalg.norm(solution['x'][-9:-11])
        return solution['info']['exitFlag'], solution['x'][-5], fDist
    
    def soc_constraints(self, tSteps, m_s, minDistance=False):
        # Creates linear, second order cone, and exponential cone constraints used 
        # in the ecos solver. These constraints are of the type, Gx <=_k h. First 
        # 2*tSteps rows in G are linear inequality constraints, the last 3*tSteps rows
        # are the exponential cones constraints, and everyting in between are SOC
        # constraints
        
        # q stores a list of the dimension of each SOC constraint, 
        # in order. Needed for ecos solver.
        
        q = [4, 3] * tSteps
        q.append(3)
        
        l = 2 * tSteps + 1          # number of linear inequalities
        e = tSteps                  # number of exponential cones
        dim_k = 7                   # number of rows in SOC constraints for one time step
        nCol  = 11 * (tSteps - 1)   # starting column for last time step
        nRow  = l + dim_k * tSteps  # starting row for final SOC constraints
        eRow  = nRow + 3            # starting row for exponetial cone constraints
        
        # the matrix G is initialized as a csc matrix, h is a vector. 
        # G_row, G_col, and G_val represent the row index, column index, and scalar
        # value, respectively, of each non-zero element in G. G and h are initially
        # set to the final mass constraint.
        
        h = np.zeros(eRow + 3*tSteps)
        
        G_row, G_col, G_val = [0], [nCol + 7], [1]
        h[0] = log(self.m_f)
        
        for t in range(tSteps):        
            k, kRow, kCol = t + 1, l + dim_k * t, 11 * t            # kRow, kCol are starting row/column for timestep
            z_0 = log(m_s - self.alpha * self.rho_2 * t * self.dt)  # z_0 is used form lower/upper thrust bounds
            
            G_row.extend([k, k])                             # upper thrust constraint
            G_col.extend([kCol + 6, kCol + 10])
            G_val.extend([-self.rho_2 * exp(-z_0), -1])
            h[k] = self.rho_2 * exp(-z_0) * (1 + z_0)
            
            G_row.extend([tSteps+k, tSteps+k])               # thrust pointing constraint
            G_col.extend([kCol + 7, kCol + 10])
            G_val.extend([1, -cos(self.theta)])
            
            G_row.extend([kRow, kRow+1, kRow+2, kRow+3])     # thrust magnitude constraint
            G_col.extend([kCol+10, kCol+7, kCol+8, kCol+9])
            G_val.extend([1, 1, 1, 1])
            
            if t != tSteps - 1:                              # glideslope constraint for every time step except
                G_row.extend(range(kRow+4, kRow+7))          # for the last time step
                G_col.extend(range(kCol, kCol+3))
                G_val.extend([1/tan(self.gamma), 1, 1])
                
                G_row.extend(range(kRow+4, kRow+7))
                G_col.extend(range(nCol, nCol+3))
                G_val.extend([-1/tan(self.gamma), -1, -1])
                    
            G_row.extend([eRow+3*t, eRow+3*t + 1])          # lower thrust bound, represented as an exponential cone
            G_col.extend([kCol + 6, kCol + 10])
            G_val.extend([-1, 1.0 / self.rho_1])
            h[eRow+3*t + 2] = 1
        
        # if golden search is being conducted, then final landing distance
        # must be less than slack variable. If no golden search, then final
        # distance must be less than distance found in golden search

        if minDistance:    
            G_row.extend([nRow, nRow + 1, nRow + 2])
            G_col.extend([11 * tSteps, nCol + 1, nCol + 2])
            G_val.extend([1, 1, 1])
        
        else:
            G_row.extend([nRow + 1, nRow + 2])
            G_col.extend([nCol + 1, nCol + 2])
            G_val.extend([1, 1]) 

            h[nRow] = self.dMax 

        G_val = [-g for g in G_val]     # g matrix has to be negated for ecos solver

        return csc_matrix((G_val, (G_row, G_col)), shape=(eRow + 3*tSteps, 11 * tSteps + 1)), h, q, l, e

    def equality_constraints(self, tSteps, x_s):   
        # Creates the equality constraints, Ax = b that will be used in 
        # the ecos solver. Defines state vector at every time step, along 
        # with initial and final conditions.
        
        n_k, nVals               = len(self.rowList), len(self.rowList) * tSteps + 16
        A_row, A_col, A_val      = np.zeros(nVals), np.zeros(nVals), np.zeros(nVals)
        A_row[0:11], A_col[0:11] = np.linspace(0, 10, 11), np.linspace(0, 10, 11)
        A_val[0:11] = np.ones(11)
        
        # Iterates through every time step (except for final time step), fills in
        # data for row, column, and value for every non-zero element in A
        
        for k in range(tSteps - 1):
            start, stop = 11 + n_k * k, 11 + n_k * (k + 1)
            
            A_row[start:stop] = [11 + a + k * 7 for a in self.rowList]
            A_col[start:stop] = [a + k * 11 for a in self.colList]
            A_val[start:stop] = self.valList
        
        # Defines final conditions. Final velocity components are all set to zero,
        # final vertical position also set to zero
        
        start, stop = 11 + n_k * (tSteps - 1), 14 + n_k * (tSteps - 1)
        nRow, nCol = 7 * (tSteps - 1) + 11, 11 * (tSteps - 1) + 3
        
        A_row[start: stop] = [nRow, nRow + 1, nRow + 2]
        A_col[start: stop] = [nCol, nCol + 1, nCol + 2]
        A_val[start: stop] = [1, 1, 1]

        A_row[stop] = nRow + 3
        A_col[stop] = nCol - 3
        A_val[stop] = 1
        
        # b represents the right hand side of the equality constraint, Ax=b. As with
        # the matrix E, bVect is repeated for every time step (except for the last time
        # step). The first 11 elements in b are the initial conditions of the system.
        
        b = np.concatenate((x_s, np.zeros(7 * tSteps)))

        for t in range(tSteps - 1):
            b[11 + t * 7: 18 + t * 7] = self.bVect
        b[11:18] = np.zeros(7)

        return csc_matrix((A_val, (A_row, A_col)), shape=(11 + 7 * tSteps, 11 * tSteps + 1)), b.astype('float64')

    