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


class OptimizeHistory(object):
    """
    
    """
    def __init__(self, param_names, objective_names):
        """
        
        :param param_names: list of str 
        :param objective_names: list of str
        """


class Individual:
    """
    
    """
    def __init__(self, x):
        """
        
        :param x: array 
        """
        self.x = x
        self.objectives = None
        self.rank = None
        self.distance = None
        self.fitness = None


class BasinhoppingGen(object):
    """
    
    """
    def __init__(self, pop_size, param_names, param_bounds, objective_names, interval=None, custom=None):
        pass


class EmooGen:
    """
    
    """
    def __init__(self, pop_size, param_names, param_bounds, objective_names, interval=None, custom=None):
        """
        
        :param pop_size: int 
        :param param_names: list of str
        :param param_bounds: list of tuple
        :param objective_names: list of str
        :param interval: int : how often to decrease the strength of mutation and crossover
        :param custom: callable, operates on a population
        """
        self.pop_size = pop_size
        self.param_names = param_names
        self.param_min = np.array([param_bound[0] for param_bound in param_bounds])
        self.param_max = np.array([param_bound[1] for param_bound in param_bounds])
        self.objective_names = objective_names
        self.custom = custom
        self.population = []
        self.config()
        self.evaluated = False

    @property
    def num_params(self):
        """
        
        :return: int 
        """
        return len(self.param_names)

    def config(self, m0=20, c0=20, p_m=0.5, delta_m=0, delta_c=0, mutate_parents=False, verbose=False):
        """
        
        :param m0: int : initial strength of mutation
        :param c0: int : initial strength of crossover
        :param p_m: float : probability of mutation  
        :param delta_m: int : decrease mutation strength every interval
        :param delta_c: int : decrease crossover strength every interval
        :param verbose: bool
        """
        self.m0 = m0
        self.m = self.m0
        self.c0 = c0
        self.c = self.c0
        self.p_m = p_m
        self.delta_m = delta_m
        self.delta_c = delta_c
        self.mutate_parents = mutate_parents
        self.verbose = verbose

    def get_random_params(self):
        """
        
        :return: array 
        """
        return np.random.uniform(self.param_min, self.param_max)

    def init_population(self):
        """
         
        """
        self.population = []
        for i in range(self.pop_size):
            params = self.random_params()
            individual = Individual(params)
            self.population.append(individual)
        self.evaluated = False
    
    def return_to_bounds(self, p):
        """
        
        :param p: array
        :return: array 
        """
        p = np.minimum(p, self.param_max)
        p = np.maximum(p, self.param_min)
        return p

    def evolve(self, maxgen=200):
        """
        Generator yields a new population. Requires that features have been evaluated and fitness assigned to current
        population.
        :param maxgen: int
        :yield: list of :class:'Individual' 
        """
        self.current_gen = 0
        while self.current_gen < maxgen:
            if self.current_gen == 0:
                self.m = self.m0
                self.c = self.c0
                self.init_population()
                if self.interval is None:
                    self.interval = maxgen
                if self.verbose:
                    print 'Starting evolutionary multiobjective optimization generator (EmooGen)\n'
                    print 'Based on Bahl A, Stemmler MB, Herz AVM, Roth A. (2012). J Neurosci Methods.\n'
                    print 'Modified by Aaron D. Milstein, Grace Ng, Ivan Soltesz (2017).'
                yield self.population
            elif not self.evaluated:
                raise Exception('EmooGen step: evolution; fitness of current population has not been evaluated.')
            else:
                if self.current_gen % self.interval == 0:
                    self.m += self.delta_m
                    self.c += self.delta_c
                    if self.verbose:
                        print 'Generation %i/%i: Decreasing strength of mutation and crossover' % \
                              (self.current_gen, maxgen)
                self.selection()
                self.crossover()
                self.mutation()
                yield self.population
            self.current_gen += 1

            # self.evaluate()
            # self.assign_fitness()
            # if (self.checkpopulation != None):
            #    self.checkpopulation(self.population)
        self.report()

    def selection(self):
        """
        In this step the mating pool is formed by selection. The population is shuffled, each individual is compared to 
        its neighbor, and the individual with high fitness score is transferred into the mating pool. This procedure is 
        repeated twice.
        """
        if not self.evaluated:
            raise Exception('EmooGen step: selection; Fitness of current population has not been evaluated.')

        mating_pool = []
        
        for k in range(1):
            population_permutation = self.population[np.random.permutation(len(self.population))]

            for i in np.arange(0, len(self.population)-1, 2):
                individual1, individual2 = population_permutation[i], population_permutation[i+1]
                if individual1.fitness < individual2.fitness:
                    mating_pool.append(individual1)
                else:
                    mating_pool.append(individual2)
        self.population = list(mating_pool)

    def crossover(self):
        """
         
        """
        children = []
        # do not add more children then original population size
        while len(children) + len(self.population) < 2 * self.pop_size:
            i, j = np.random.choice(range(len(self.population)), 2)
            parent1 = self.population[i]
            parent2 = self.population[j]
            child1_params = np.empty(self.num_params)
            child2_params = np.empty(self.num_params)
            for i in range(self.num_params):
                u_i = np.random.random()
                if u_i <= 0.5:
                    beta_q_i = pow(2. * u_i, 1. / (self.c + 1))
                else:
                    beta_q_i = pow(1. / (2. * (1. - u_i)), 1. / (self.c + 1))
                child1_params[i] = 0.5 * ((1. + beta_q_i) * parent1.p[i] + (1. - beta_q_i) * parent2.p[i])
                child2_params[i] = 0.5 * ((1. - beta_q_i) * parent1.p[i] + (1 + beta_q_i) * parent2.p[i])
            child1 = Individual(self.return_to_bounds(child1_params))
            child2 = Individual(self.return_to_bounds(child2_params))
            children.append(child1)
            children.append(child2)
        self.population.extend(children)

    def mutation(self):
        """
        polynomial mutation (Deb, 2001)
        """
        for k in range(len(self.population)):
            individual = self.population[k]
            if self.mutate_parents or individual.fitness is None:
                individual.fitness = None
                for i in range(self.num_params):
                    # each gene only mutates with a certain probability
                    if np.random.random() < self.p_m:
                        r_i = np.random.random()
                        if r_i < 0.5:
                            delta_i = pow(2. * r_i, 1. / (self.m + 1)) - 1.
                        else:
                            delta_i = 1. - pow(2. * (1. - r_i), 1. / (self.m + 1))
                        individual.p[i] += delta_i
                individual.p = self.return_to_bounds(individual.p)

    def evaluate(self):
        # only evaluate up to pop_size, as that number of processes must be pre-allocated
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
