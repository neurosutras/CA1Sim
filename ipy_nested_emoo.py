"""
(C) Armin Bahl 16.01.2009, UCL, London, UK
modified on ACCN 2011 Bedlewo, Poland 20.08.2011 
further modification: 23.04.2012 (Munich)

If you use this algorithm for your research please cite:

Bahl A, Stemmler MB, Herz AVM, Roth A. (2012). Automated
optimization of a reduced layer 5 pyramidal cell model based on
experimental data. J Neurosci Methods; in press
"""
import numpy as np
import time

try:
    from mpi4py import MPI
    mpi4py_loaded = True
except:
    mpi4py_loaded = False


class NestedEmoo:
     
    def __init__(self, N, C, variables, objectives, infos=[]):
        self.version = 1.0
        # standard population parameters
        self.size = N              # size of population, must be even and a multiple of processors - 1
        self.capacity = C      # when new children are born, we have this amount of individuals
        
        # number of processors
        self.variables = variables
        self.obj = len(objectives)                # number of objectives
        self.infos = len(infos)
        
        self.objectives_names = objectives
        self.infos_names = infos
        
        self.para = len(self.variables)
                
        self.no_properties = np.ones(3)*(-1.0)
        self.no_objectives = np.ones(self.obj+self.infos)*(-1)
        
        self.columns = dict({})
        self.column_names = []
        
        
        self.objpos = self.para
        self.infopos = self.objpos + self.obj
        self.rankpos = self.infopos + self.infos
        self.distpos = self.rankpos + 1
        self.fitnesspos = self.distpos + 1
            
        i = 0
        for variable in variables:
            self.column_names.append(variable[0])
            self.columns[variable[0]] = i
            i+=1
            
        for objective in objectives:
            self.column_names.append(objective)
            self.columns[objective] = i
            i+=1
        
        for info in infos:
            self.column_names.append(info)
            self.columns[info] = i
            i+=1
            
        self.column_names.append('emoo_rank')
        self.columns['emoo_rank'] = self.rankpos
        self.column_names.append('emoo_dist')
        self.columns['emoo_dist'] = self.distpos
        self.column_names.append('emoo_fitness')
        self.columns['emoo_fitness'] = self.fitnesspos 
        
        self.checkfullpopulation = None
        self.checkpopulation = None
        self.setuped = False
        
        if mpi4py_loaded == True:
            self.comm = MPI.COMM_WORLD
            self.master_mode = self.comm.rank == 0
            self.mpi = self.comm.size > 1
        else:
            self.master_mode = True
            self.mpi = False
        
        
    def setup(self, eta_m_0=20, eta_c_0=20, p_m=0.5, finishgen=-1, d_eta_m=0, d_eta_c=0, mutate_parents=False):
        self.eta_m_0 = eta_m_0
        self.eta_c_0 = eta_c_0
        self.p_m = p_m
        self.finishgen = finishgen
        self.d_eta_m = d_eta_m
        self.d_eta_c = d_eta_c
        self.mutate_parents = mutate_parents
        
        self.setuped = True
        
    def normit(self, p):
        p_norm = np.zeros(len(p), dtype=float)
        
        for i in range(len(p)):
            p_norm[i] = (p[i]-self.variables[i][1])/(self.variables[i][2] - self.variables[i][1])
                   
        return p_norm

    def unnormit(self, p_norm):
        p = np.zeros(len(p_norm), dtype=float)
    
        for i in range(len(p_norm)):
            p[i] = p_norm[i]*(self.variables[i][2] - self.variables[i][1]) + self.variables[i][1]
                   
        return p
    
    def getpopulation_unnormed(self):
        unnormed_population = []
        for individual in self.population:
            individual_unnormed = individual.copy()
            individual_unnormed[:self.para] = self.unnormit(individual[:self.para])
            unnormed_population.append(individual_unnormed)

        return np.array(unnormed_population)

    def initpopulation(self):
        
        init_parameters = np.random.rand(self.size, self.para)
        init_properties = np.ones((self.size, self.obj+self.infos+3))*(-1.0)

        self.population = np.c_[init_parameters, init_properties]       
    
    def evolution(self, generations):
        if self.setuped == False:
            print "Please run setup"
            return
        
        if(self.master_mode == True):
            
            self.eta_c = self.eta_c_0
            self.eta_m = self.eta_m_0
            
            self.initpopulation()
            
            print "This is Evolutionary Multiobjective Optimization (Emoo), Version %.1f."%self.version
            print "\n www.g-node.org/emoo" 
            print "\n\n\tIf you use this algorithm for your research please cite:"
            print "\n\t\tBahl A, Stemmler MB, Herz AVM, Roth A. (2012). Automated optimization of a reduced"
            print "\t\tlayer 5 pyramidal cell model based on experimental data. J Neurosci Methods.\n"
            print "\n\t ... And please provide us with feedback. Thanks!"
            if self.mpi:
                print "\nRunning Emoo on %d processors"%self.comm.size
                print " ... let the nodes startup. Starting Optimization in 5 seconds..."
                time.sleep(5) # let all the slaves load
            
        
            print "Starting Evolution..."
            
            self.evaluate()

            self.assign_fitness()
            
            if(self.checkpopulation != None):
                self.checkpopulation(self.getpopulation_unnormed(), self.columns, 0)
                    
            for gen in range(1, generations):
                
                # Change the Crossover and Mutation Parameters
                if gen > self.finishgen and self.finishgen != -1:
                    self.eta_c += self.d_eta_c
                    self.eta_m += self.d_eta_m
                      
                self.selection()
                self.crossover()
                self.mutation()
                
                self.evaluate()

                self.assign_fitness()
                
                if self.checkfullpopulation != None:
                    self.checkfullpopulation(self.getpopulation_unnormed(),self.columns, gen)
                    
                self.new_generation()
                
                if(self.checkpopulation != None):
                    self.checkpopulation(self.getpopulation_unnormed(), self.columns, gen)
        
            # tell the slaves (if any) to terminate
            if self.mpi == True:
                for i in range(1, self.comm.size):
                    self.comm.send(None, dest=i)
                
                time.sleep(5) # let all the slaves finish
            
            print "Evolution done!!!"

        else:      
            self.evaluate_slave()
                
                
    def selection(self):
        """
        In this step the mating pool is formed by selection
        The population is shuffelded and then each individal is compared with the next and only
        the better will be tranfered into the mating pool
        then the population is shuffelded again and the same happens again
        """
        
        # the population has the size N now
        # and all fitnesses are assigned!
        
        mating_pool = []
        
        for k in [0,1]:
            population_permutation = self.population[np.random.permutation(len(self.population))]
            # -1 because in the cases off odd population size!
            for i in np.arange(0, len(self.population)-1, 2):
                fitness1 = population_permutation[i][-1]
                fitness2 = population_permutation[i+1][-1]
                
                if(fitness1 < fitness2):
                    mating_pool.append(population_permutation[i])
                else:
                    mating_pool.append(population_permutation[i+1])
        
        # now we have a mating pool
        
        # this is our new population
        self.population = np.array(mating_pool)         
            
    def crossover(self):
        
        children = []
        
        while(len(children) + len(self.population) < self.capacity):
            
            # choose two random parents
            p = int(np.random.random()*len(self.population))
            q = int(np.random.random()*len(self.population))
            
            parent1 = self.population[p][:self.para]
            parent2 = self.population[q][:self.para]
            
            parameters1 = np.empty(self.para)
            parameters2 = np.empty(self.para)
                
            # determine the crossover parameters
            for i in range(self.para):
                u_i = np.random.random()
            
                if u_i <= 0.5:
                    beta_q_i = pow(2.*u_i, 1./(self.eta_c+1))
                else:
                    beta_q_i = pow(1./(2*(1-u_i)), 1./(self.eta_c+1))
            
                parameters1[i]  = 0.5*((1+beta_q_i)*parent1[i] + (1-beta_q_i)*parent2[i])
                parameters2[i]  = 0.5*((1-beta_q_i)*parent1[i] + (1+beta_q_i)*parent2[i])
            
                # did we leave the boundary?
                if(parameters1[i] > 1):
                    parameters1[i] = 1
                
                if(parameters1[i] < 0):
                    parameters1[i] = 0
                
                if(parameters2[i] > 1):
                    parameters2[i] = 1
                
                if(parameters2[i] < 0):
                    parameters2[i] = 0
                
            offspring1 = np.r_[parameters1, self.no_objectives, self.no_properties]
            offspring2 = np.r_[parameters2, self.no_objectives, self.no_properties]

            children.append(offspring1)
            children.append(offspring2)               

        children = np.array(children)
        self.population = np.r_[self.population, children]

         
    def mutation(self):
        
        # polynomial mutation (Deb, 124)
        for k in range(len(self.population)):
            
            individual = self.population[k]
            
            if not self.mutate_parents and individual[self.fitnesspos] != -1:
                continue # this is a parent, do not mutate it
            
            for i in range(self.para):
            
                # each gene only mutates with a certain probability
                m = np.random.random()
                
                if(m < self.p_m):
                    r_i = np.random.random()
                
                    if r_i < 0.5:
                        delta_i = pow(2*r_i, 1./(self.eta_m+1)) - 1
                    else:
                        delta_i = 1-pow(2*(1-r_i), 1./(self.eta_m+1))
                        
                    individual[i] += delta_i
                    
                    # did we leave the boundary?
                    if(individual[i] > 1):
                        individual[i] = 1
                    
                    if(individual[i] < 0):
                        individual[i] = 0
            
            individual[self.para:] = np.r_[self.no_objectives, self.no_properties]

    def evaluate(self):
        
        new_population = []
        
        # is the master alone?
        if(self.mpi == False):

            for individual in self.population:
                
                # only evaluate those that are really new!
                if individual[self.fitnesspos] == -1:
                    
                    parameters = individual[:self.para]
                   
                    objectives_error = self.evaluate_individual(parameters)
                   
                    if(objectives_error != None):
                        new_population.append(np.r_[parameters, objectives_error, self.no_properties])
                else:
                    new_population.append(individual)
        else:
            # distribute the individuals among the slaves
            i = 0
            for individual in self.population:
                if individual[self.fitnesspos] == -1:
                    parameters = individual[:self.para]
                
                    dest = i%(self.comm.size-1) + 1
                    self.comm.send(parameters, dest=dest)
                    i += 1
                else:
                    new_population.append(individual)
            
            # the master does also one
            # TODO
            
            # Receive the results from the slaves
            for i in range(i):
                result = self.comm.recv(source=MPI.ANY_SOURCE)
                
                if result != None:
                    new_population.append(np.r_[result[0], result[1], self.no_properties])
        
        self.population = np.array(new_population)
    
    def evaluate_individual(self, parameters):
        
        parameters_unnormed = self.unnormit(parameters)
                    
        # make a dictionary with the unormed parameters and send them to the evaluation function
        dict_parameters_normed = dict({})
        for i in range(len(self.variables)):
            dict_parameters_normed[self.variables[i][0]] = parameters_unnormed[i]
        
        dict_results = self.get_objectives_error(dict_parameters_normed)
        
        list_results = []
        for objective_name in self.objectives_names:
            list_results.append(dict_results[objective_name])
            
        for info_name in self.infos_names:
            list_results.append(dict_results[info_name])
        
        return np.array(list_results)
        
    def evaluate_slave(self):
        
        # We wait for parameters
        # we do not see the whole population!
        
        while(True):
            parameters = self.comm.recv(source=0) # wait....
            
            # Does the master want the slave to shutdown?
            if(parameters == None):
                # Slave finishing...
                break
            
            objectives_error = self.evaluate_individual(parameters)
            
            #objectives_error = self.get_objectives_error(self.unnormit(parameters))
            if(objectives_error == None):
                self.comm.send(None, dest=0)
            else: 
                self.comm.send([parameters, objectives_error], dest=0)
    
    def assign_fitness(self):           
        """
        are we in a multiobjective regime, then the selection of the best individual is not trival
        and must be based on dominance, thus we determine all non dominated fronts and only use the best
        to transfer into the new generation
        """
        if(self.obj > 1):
            self.assign_rank()

            new_population = np.array([])
            
            maxrank = self.population[:,self.rankpos].max()
            
            for rank in range(0, int(maxrank)+1):
                
                new_front = self.population[np.where(self.population[:,self.rankpos] == rank)]
                
                new_sorted_front = self.crowding_distance_sort(new_front)
                
                if(len(new_population) == 0):
                    new_population = new_sorted_front
                else:
                    new_population = np.r_[new_population, new_sorted_front]
                
            self.population = new_population
                         
        else:
            # simple sort the objective value
            ind = np.argsort(self.population[:,self.objpos])
            self.population = self.population[ind]
        
        # now set the fitness, indiviauls are sorted, thus fitnes is easy to set
        fitness = range(0, len(self.population[:,0]))
        self.population[:,-1] = fitness   
                    
    
    def new_generation(self):
        # the worst are at the end, let them die, if there are too many
        if(len(self.population) > self.size):
            self.population = self.population[:self.size]
         
    def dominates(self, p, q):
        
        objectives_error1 = self.population[p][self.objpos:self.objpos+self.obj]
        objectives_error2 = self.population[q][self.objpos:self.objpos+self.obj]
        
        diff12 = objectives_error1 - objectives_error2
        
        # is individdum equal or better then individdum two?
        # and at least in one objective better
        # then it dominates individuum2
        # if not it does not dominate two (which does not mean that 2 may not dominate 1)
        return ( ((diff12<= 0).all()) and ((diff12 < 0).any()) )

    
    def assign_rank(self):
            
        F = dict()

        P = self.population
        
        S = dict()
        n = dict()
        F[0] = []
        
        # determine how many solutions are dominated or dominate
        for p in range(len(P)):
            
            S[p] = []       # this is the list of solutions dominated by p
            n[p] = 0        # how many solutions are dominating p
            
            for q in range(len(P)):
                
                if self.dominates(p, q):
                    S[p].append(q)      # add q to the list of solutions dominated by p
                elif self.dominates(q, p):
                    n[p] += 1           # q dominates p, thus increase number of solutions that dominate p
                
            
            if n[p] == 0:       # no other solution dominates p
                
                # this is the rank column
                P[p][self.rankpos] = 0
                
                F[0].append(p)  # add p to the list of the first front
            
        # find the other non dominated fronts
        i = 0
        while len(F[i]) > 0:
            Q = []              # this will be the next front
            
            # take the elements from the last front
            for p in F[i]:
                
                # and take the elements that are dominated by p
                for q in S[p]:
                    # decrease domination number of all elements that are dominated by p
                    n[q] -= 1
                    # if the new domination number is zero, than we have found the next front       
                    if n[q] == 0:
                        
                        P[q][self.rankpos] = i + 1
                        Q.append(q)
            
            i += 1
            F[i] = Q    # this is the next front
    
    
    def crowding_distance_sort(self, front):
        
        sorted_front = front.copy()
        
        l = len(sorted_front[:,0])
        
        sorted_front[:,self.distpos] = np.zeros_like(sorted_front[:,0])
        
        for m in range(self.obj):
            ind = np.argsort(sorted_front[:,self.objpos + m])
            sorted_front = sorted_front[ind]
            
            # definitely keep the borders
            sorted_front[0, self.distpos] += 1000000000000000.
            sorted_front[-1, self.distpos] += 1000000000000000.

            fm_min = sorted_front[0, self.objpos + m]
            fm_max = sorted_front[-1, self.objpos + m]
            
            if fm_min != fm_max:
                for i in range(1, l - 1):
                    sorted_front[i, self.distpos] += (sorted_front[i+1, self.objpos + m] - sorted_front[i-1, self.objpos + m])/(fm_max - fm_min)

        ind = np.argsort(sorted_front[:,self.distpos])
        sorted_front = sorted_front[ind]
        sorted_front = sorted_front[-1 - np.arange(len(sorted_front))]
                                                         
        return sorted_front

        
