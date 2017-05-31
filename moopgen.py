__author__ = 'Aaron D. Milstein'

import numpy as np
import scipy.optimize
import collections
from scipy._lib._util import check_random_state
from copy import deepcopy

"""
Here we have used scipy.optimize.basinhopping and emoo as inspiration in designing a parallel computing-compatible
framework for multi-objective parameter optimization with bounds. Rather than processing and storing a single parameter 
array at a time, these classes contain iterators and evaluation methods to process many parameter arrays in parallel,
and to store a complete history for later inspection .
"""

#__all__ = ['BGen']


class Individual(object):
    """

    """
    def __init__(self, x):
        """

        :param x: array
        """
        self.x = np.array(x)
        self.objectives = None
        self.energy = None
        self.rank = None
        self.distance = None
        self.fitness = None
        self.survivor = False


class PopulationStorage(object):
    """
    Class used to store populations of parameters and objectives during optimization.
    """
    def __init__(self, param_names, objective_names):
        """

        :param param_names:
        :param objective_names:
        """
        self.param_names = param_names
        self.objective_names = objective_names
        self.history = []

    def append(self, population):
        """

        :param population: list of :class:'Individual'
        """
        self.history.append(deepcopy(population))

    def get_best(self, n=1, generation=None, evaluate=None):
        """
        If 'last' generation is specified, and rankings have not already been stored, compute new rankings.
        If generations is specified as an integer q, compute new rankings for the last q generations.
        If 'all' generations is specified, collapse across all generations, exclude copies of Individuals that survived
        across generations, and compute new global rankings.
        Return the n best.
        :param n: int
        :param generation: str or int
        :param evaluate: callable
        :return: list of :class:'Individual'
        """
        if generation is None:
            generation = 'all'
        elif generation not in ['all', 'last'] and type(generation) != int:
            print 'PopulationStorage: Defaulting to get_best across all generations.'
            generation = 'all'
        if evaluate is None:
            evaluate = evaluate_basinhopping
        elif not isinstance(evaluate, collections.Callable):
            raise TypeError("PopulationStorage: evaluate must be callable.")
        if generation == 'last':
            group = [individual for individual in self.history[-1] if individual.rank is not None]
            if not group:
                evaluate(self.history[-1])
                group = [individual for individual in self.history[-1] if individual.rank is not None]
        elif generation == 'all':
            group = [deepcopy(individual) for population in self.history for individual in population
                     if not individual.survivor]
            evaluate(group)
        else:
            group = [deepcopy(individual) for population in self.history[-generation:] for individual in population]
            evaluate(group)
        indexes = range(len(group))
        rank = [individual.rank for individual in group]
        indexes.sort(key=rank.__getitem__)
        group = map(group.__getitem__, indexes)

        return group[:n]


class BoundedStep(object):
    """
    Step-taking method for use with BGen. Steps each parameter within specified bounds. Explores the range in log10
    space when the range is greater than 2 orders of magnitude. If bounds are not provided for some parameters, the
    default is 0.1 and 10. * x0.
    """
    def __init__(self, x0, bounds=None, stepsize=0.5, random=None):
        """

        :param x0: array
        :param bounds: list of tuple
        :param stepsize: float
        :param random: int or :class:'np.random.RandomState'
        """
        self.stepsize = stepsize
        if random is None:
            self.random = np.random
        else:
            self.random = random
        if bounds is None:
            xmin = None
            xmax = None
        else:
            xmin = [bound[0] for bound in bounds]
            xmax = [bound[1] for bound in bounds]
        if xmin is None:
            xmin = [None for i in range(len(x0))]
        if xmax is None:
            xmax = [None for i in range(len(x0))]
        for i in range(len(x0)):
            if xmin[i] is None:
                if x0[i] > 0.:
                    xmin[i] = 0.1 * x0[i]
                elif x0[i] == 0.:
                    xmin[i] = -1.
                else:
                    xmin[i] = 10. * x0[i]
            if xmax[i] is None:
                if x0[i] > 0.:
                    xmax[i] = 10. * x0[i]
                elif x0[i] == 0.:
                    xmax[i] = 1.
                else:
                    xmax[i] = 0.1 * x0[i]
        self.x0 = x0
        self.x_range = np.subtract(xmax, xmin)
        self.order_mag = np.ones_like(x0)
        if not np.any(np.array(xmin) == 0.):
            self.order_mag = np.abs(np.log10(np.abs(np.divide(xmax, xmin))))
        else:
            for i in range(len(x0)):
                if xmin[i] == 0.:
                    self.order_mag[i] = int(xmax[i] / 10)
                else:
                    self.order_mag[i] = abs(np.log10(abs(xmax[i] / xmin[i])))
        self.log10_range = np.log10(np.add(1., self.x_range))
        self.x_offset = np.subtract(1., xmin)

    def __call__(self, current_x):
        """

        :param current_x: array
        :return: array
        """
        x = np.add(current_x, self.x_offset)
        x = np.maximum(x, 1.)
        x = np.minimum(x, np.add(1., self.x_range))
        for i in range(len(x)):
            if self.order_mag[i] >= 2.:
                x[i] = self.log10_step(i, x[i])
            else:
                x[i] = self.linear_step(i, x[i])
        new_x = np.subtract(x, self.x_offset)
        return new_x

    def linear_step(self, i, xi):
        """

        :param i: int
        :param xi: float
        :return: float
        """
        step = self.stepsize * self.x_range[i] / 2.
        new_xi = self.random.uniform(max(1., xi-step), min(xi+step, 1.+self.x_range[i]))
        return new_xi

    def log10_step(self, i, xi):
        """

        :param i: int
        :param xi: float
        :return: float
        """
        step = self.stepsize * self.log10_range[i] / 2.
        xi = np.log10(xi)
        new_xi = self.random.uniform(max(0., xi-step), min(xi+step, self.log10_range[i]))
        new_xi = np.power(10., new_xi)
        return new_xi


def sort_by_crowding_distance(population):
    """
    Modifies in place the distance attribute of each Individual in the population. Returns the sorted population.
    :param population: list of :class:'Individual'
    :return: list of :class:'Individual'
    """
    pop_size = len(population)
    num_objectives = [len(individual.objectives) for individual in population if individual.objectives is not None]
    if len(num_objectives) < pop_size:
        raise Exception('sort_by_crowding_distance: objectives have not been stored for all Individuals in population')
    num_objectives = max(num_objectives)
    for individual in population:
        individual.distance = 0
    for m in range(num_objectives):
        indexes = range(pop_size)
        objective_vals = [individual.objectives[m] for individual in population]
        indexes.sort(key=objective_vals.__getitem__)
        population = map(population.__getitem__, indexes)

        # keep the borders
        population[0].distance += 1.e15
        population[-1].distance += 1.e15

        objective_min = population[0].objectives[m]
        objective_max = population[-1].objectives[m]

        if objective_min != objective_max:
            for i in range(1, pop_size - 1):
                population[i].distance += (population[i + 1].objectives[m] - population[i - 1].objectives[m]) / \
                                     (objective_max - objective_min)
    indexes = range(pop_size)
    distances = [individual.distance for individual in population]
    indexes.sort(key=distances.__getitem__)
    indexes.reverse()
    population = map(population.__getitem__, indexes)
    return population


def assign_rank_by_energy(population):
    """
    Modifies in place the rank attributes of each Individual in the population. Sorts the population by total energy
    (sum of all objectives).
    :param population: list of :class:'Individual'
    """
    pop_size = len(population)
    indexes = range(len(population))
    energy_vals = [individual.energy for individual in population if individual.energy is not None]
    if len(energy_vals) < pop_size:
        raise Exception('assign_rank_by_energy: energy has not been stored for all Individuals in population')
    indexes.sort(key=energy_vals.__getitem__)
    population = map(population.__getitem__, indexes)
    for rank, individual in enumerate(population):
        individual.rank = rank


def assign_rank_by_crowding_distance(population):
    """
    Modifies in place the distance and rank attributes of each Individual in the population. This is appropriate for
    early generations of evolutionary optimization, and helps to preserve diversity of solutions. However, once all
    members of the population have converged to a single fitness value, naive ranking by crowding distance can favor
    unique solutions over lower energy solutions. In this case, rank is assigned by total energy.
    :param population: list of :class:'Individual'
    """
    pop_size = len(population)
    fitness_vals = [individual.fitness for individual in population]
    if len(fitness_vals) < pop_size:
        raise Exception('assign_rank_by_crowding_distance: fitness has not been stored for all Individuals in '
                        'population')
    max_fitness = max(fitness_vals)
    if max_fitness > 0:
        new_population = []
        for fitness in range(max_fitness + 1):
            new_front = [individual for individual in population if individual.fitness == fitness]
            new_sorted_front = sort_by_crowding_distance(new_front)
            new_population.extend(new_sorted_front)
        # now that population is sorted, assign rank to Individuals
        for rank, individual in enumerate(new_population):
            individual.rank = rank
    else:
        assign_rank_by_energy(population)


def assign_fitness_by_dominance(population, disp=False):
    """
    Modifies in place the fitness attribute of each Individual in the population.
    :param population: list of :class:'Individual'
    :param disp: bool
    """
    def dominates(p, q):
        """
        Individual p dominates Individual q if each of its objective values is equal or better, and at least one of
        its objective values is better.
        :param p: :class:'Individual'
        :param q: :class:'Individual'
        :return: bool
        """
        diff12 = np.subtract(p.objectives, q.objectives)
        return ((diff12 <= 0.).all()) and ((diff12 < 0.).any())

    pop_size = len(population)
    num_objectives = [len(individual.objectives) for individual in population if individual.objectives is not None]
    if len(num_objectives) < pop_size:
        raise Exception('assign_fitness_by_dominance: objectives have not been stored for all Individuals in '
                        'population')
    num_objectives = max(num_objectives)
    if num_objectives > 1:
        F = {0: []}  # first front of dominant Individuals
        S = dict()
        n = dict()

        for p in range(len(population)):
            S[p] = []  # list of Individuals that p dominates
            n[p] = 0  # number of Individuals that dominate p

            for q in range(len(population)):
                if dominates(population[p], population[q]):
                    S[p].append(q)
                elif dominates(population[q], population[p]):
                    n[p] += 1

            if n[p] == 0:
                population[p].fitness = 0  # fitness 0 indicates first dominant front
                F[0].append(p)

        # excluding the Individuals that dominated the previous front, find the next front
        i = 0
        while len(F[i]) > 0:
            F[i+1] = []  # next front
            # take the elements from the previous front
            for p in F[i]:
                # take the elements that p dominates
                for q in S[p]:
                    # decrease domination value of all Individuals that p dominates
                    n[q] -= 1
                    if n[q] == 0:
                        population[q].fitness = i + 1  # assign fitness of current front
                        F[i+1].append(q)
            i += 1
    else:
        for individual in population:
            individual.fitness == 0
    if disp:
        print F


def evaluate_basinhopping(population, disp=False):
    """
    Modifies in place the fitness and rank attributes of each Individual in the population.
    :param population: list of :class:'Individual'
    :param disp: bool
    :return:
    """
    assign_fitness_by_dominance(population)
    assign_rank_by_energy(population)


def evaluate_random(population, disp):
    """
    Modifies in place the rank attribute of each Individual in the population.
    :param population: list of :class:'Individual'
    """
    rank_vals = range(len(population))
    np.random.shuffle(rank_vals)
    for i, individual in enumerate(population):
        rank = rank_vals[i]
        individual.rank = rank
        if disp:
            print 'Individual %i: rank %i, x: %s' % (i, rank, individual.x)


class BGen(object):
    """
    The class is inspired by scipy.optimize.basinhopping. It provides a generator interface to produce a list of
    parameter arrays for parallel evaluation. Features adaptive pruning and reduction of step_size every iteration.
    """
    def __init__(self, x0, param_names, objective_names, pop_size, bounds=None, take_step=None, evaluate=None,
                 seed=None, max_gens=100, adaptive_step_interval=20, adaptive_step_factor=0.9, ngen_success=None,
                 survival_rate=0.1, disp=False):
        """

        :param x0: array
        :param param_names: list of str
        :param objective_names: list of str
        :param pop_size: int
        :param bounds: list of tuple
        :param take_step: callable
        :param evaluate: callable
        :param seed: int or :class:'np.random.RandomState'
        :param max_gens: int
        :param adaptive_step_interval: int
        :param adaptive_step_factor: float in [0., 1.]
        :param ngen_success: int
        :param survival_rate: float in [0., 1.]
        :param disp: bool
        """
        self.x0 = np.array(x0)
        self.storage = PopulationStorage(param_names, objective_names)
        self.pop_size = pop_size
        if evaluate is None:
            self._evaluate = evaluate_basinhopping
        elif isinstance(evaluate, collections.Callable):
            self._evaluate = evaluate
        else:
            raise TypeError("BGen: evaluate must be callable.")
        self.random = check_random_state(seed)
        if take_step is None:
            self.take_step = BoundedStep(x0, stepsize=0.5, bounds=bounds, random=self.random)
        elif isinstance(take_step, collections.Callable):
            self.take_step = take_step
        else:
            raise TypeError("BGen: take_step must be callable.")
        self.max_gens = max_gens
        self.adaptive_step_interval = adaptive_step_interval
        self.adaptive_step_factor = adaptive_step_factor
        if ngen_success is None:
            self.ngen_succcess = self.max_gens
        else:
            self.ngen_succcess = ngen_success
        self.num_survivors = int(pop_size * survival_rate)
        self.disp = disp

        self.objectives_stored = False
        self.evaluated = False
        self.population = []
        self.survivors = []

    def __call__(self):
        """
        A generator that yields a population of size pop_size as a list of individuals.
        :yields: list of :class:'Individual'
        """
        self.num_gen = 0
        while self.num_gen < self.max_gens:
            if self.num_gen == 0:
                self.init_population()
                self.evaluated = False
            else:
                if not self.objectives_stored:
                    raise Exception('BGen: Gen %i, objectives have not been stored for all Individuals in '
                                    'population' % self.num_gen)
                if self.num_gen % self.adaptive_step_interval == 0:
                    self.evaluate()
                    new_step_size = self.take_step.stepsize * self.adaptive_step_factor
                    if self.disp:
                        print 'BGen: Gen %i, previous step_size: %.2f, new step_size: %.2f' % \
                              (self.num_gen, self.take_step.stepsize, new_step_size)
                    self.take_step.stepsize = new_step_size
                    self.step_survivors()
                else:
                    self.step_population()
                self.evaluated = False
            self.objectives_stored = False
            if self.disp:
                print 'BGen: Gen %i, yielding parameters for population size %i' % (self.num_gen, self.pop_size)
            yield [individual.x for individual in self.population]
            self.num_gen += 1

    def set_objectives(self, objectives):
        """
        Expects a list of objective arrays to be in the same order as the list of parameter arrays yielded from the
        current generation.
        :param objectives: list of array
        """
        for i, objective_dict in enumerate(objectives):
            if type(objective_dict) != dict:
                raise TypeError('BGen.set_objectives: objectives must be a list of dict')
            this_objectives = np.array([objective_dict[key] for key in self.storage.objective_names])
            self.population[i].objectives = this_objectives
            self.population[i].energy = np.sum(this_objectives)
        if self.disp:
            print 'BGen: Gen %i, storing objectives for population size %i' % (self.num_gen, self.pop_size)
        self.objectives_stored = True

    def evaluate(self):
        """
        Assign fitness
        """
        self._evaluate(self.survivors + self.population)
        self.storage.append(self.survivors + self.population)
        self.evaluated = True

    def init_population(self):
        """

        """
        self.population = [Individual(self.take_step(self.x0)) for i in range(self.pop_size)]

    def step_survivors(self):
        """
        Consider the highest ranked Individuals of the current generation to be the survivors. Seed the next
        generation with steps taken from the surviving set of parameters.
        """
        if not self.evaluated:
            raise Exception('BGen: Current generation has not yet been evaluated.')
        self.survivors = self.storage.get_best(self.num_survivors, generation=self.adaptive_step_interval)
        for individual in self.survivors:
            individual.survivor = True
        new_population = []
        for i in range(self.pop_size):
            individual = Individual(self.take_step(self.survivors[i % self.num_survivors].x))
            new_population.append(individual)
        self.population = new_population

    def step_population(self):
        """

        """
        new_population = [Individual(self.take_step(individual.x)) for individual in self.population]
        self.population = new_population
