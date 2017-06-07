__author__ = 'Aaron D. Milstein'

from function_lib import *
from matplotlib.pyplot import cm
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
        self.features = None
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
    def __init__(self, param_names=None, feature_names=None, objective_names=None, file_path=None):
        """

        :param param_names: list of str
        :param feature_names: list of str
        :param objective_names: list of str
        :param file_path: str (path)
        """
        if param_names is None:
            self.param_names = []
        else:
            self.param_names = param_names
        if feature_names is None:
            self.feature_names = []
        else:
            self.feature_names = feature_names
        if objective_names is None:
            self.objective_names = []
        else:
            self.objective_names = objective_names
        self.history = []  # a list of populations, of length max_gens
        self.survivors = []  # a list of populations (some may be empty), of length max_gens
        if file_path is not None:
            self.load(file_path)

    def append(self, population, survivors=None):
        """

        :param population: list of :class:'Individual'
        :param survivors: list of :class:'Individual'
        """
        if survivors is None:
            survivors = []
        self.survivors.append(deepcopy(survivors))
        self.history.append(deepcopy(population))

    def get_best(self, n=1, generations=None, evaluate=None, modify=False):
        """
        If 'last' generation is specified, and rankings have not already been stored, compute new rankings.
        If generations is specified as an integer q, compute new rankings for the last q generations, including the set
        of survivors produced closest to, but before the qth generation.
        If 'all' generations is specified, collapse across all generations, exclude copies of Individuals that survived
        across generations, and compute new global rankings.
        Return the n best.
        If modify is True, allow changes to the rankings of Individuals stored in history, otherwise operate on and
        discard copies.
        :param n: int
        :param generations: str or int
        :param evaluate: callable
        :param modify: bool
        :return: list of :class:'Individual'
        """
        if generations is None:
            generations = 'all'
            print 'PopulationStorage: Defaulting to get_best across all generations.'
        elif generations not in ['all', 'last'] and type(generations) != int:
            generations = 'all'
            print 'PopulationStorage: Defaulting to get_best across all generations.'
        if evaluate is None:
            evaluate = evaluate_bgen
        elif not isinstance(evaluate, collections.Callable):
            raise TypeError("PopulationStorage: evaluate must be callable.")
        if generations == 'last':
            recent_survivors = self.survivors[-1]  # may be empty
            group = [individual for individual in self.history[-1] + recent_survivors if individual.rank is not None]
            if len(group) < len(self.history[-1] + recent_survivors):
                if modify:
                    group = [individual for individual in self.history[-1] + recent_survivors]
                else:
                    group = [deepcopy(individual) for individual in self.history[-1] + recent_survivors]
                evaluate(group)
        elif generations == 'all':
            if modify:
                group = [individual for population in self.history for individual in population]
            else:
                group = [deepcopy(individual) for population in self.history for individual in population]
            evaluate(group)
        else:
            generations = min(len(self.history), generations)
            if modify:
                group = [individual for population in self.history[-generations:] for individual in population]
                group.extend([individual for individual in self.survivors[-generations]])
            else:
                group = [deepcopy(individual) for population in self.history[-generations:] for individual in population]
                group.extend([deepcopy(individual) for individual in self.survivors[-generations]])
            evaluate(group)
        indexes = range(len(group))
        rank = [individual.rank for individual in group]
        indexes.sort(key=rank.__getitem__)
        group = map(group.__getitem__, indexes)
        return group[:n]

    def plot(self):
        """

        """
        colors = list(cm.rainbow(np.linspace(0, 1, len(self.history))))
        for this_attr in ['fitness', 'energy', 'distance', 'survivor']:
            plt.figure()
            for j, population in enumerate(self.history):
                plt.scatter([indiv.rank for indiv in population], [getattr(indiv, this_attr) for indiv in population],
                            c=colors[j], alpha=0.1)
                plt.scatter([indiv.rank for indiv in self.survivors[j]],
                            [getattr(indiv, this_attr) for indiv in self.survivors[j]], c=colors[j], alpha=0.5)
            plt.title(this_attr)
        for i, param_name in enumerate(self.param_names):
            this_attr = 'x'
            plt.figure()
            for j, population in enumerate(self.history):
                plt.scatter([indiv.rank for indiv in population],
                            [getattr(indiv, this_attr)[i] for indiv in population],
                            c=colors[j], alpha=0.1)
                plt.scatter([indiv.rank for indiv in self.survivors[j]],
                            [getattr(indiv, this_attr)[i] for indiv in self.survivors[j]], c=colors[j], alpha=0.5)
            plt.title(param_name)
        for i, objective_name in enumerate(self.objective_names):
            this_attr = 'objectives'
            plt.figure()
            for j, population in enumerate(self.history):
                plt.scatter([indiv.rank for indiv in population],
                            [getattr(indiv, this_attr)[i] for indiv in population],
                            c=colors[j], alpha=0.1)
                plt.scatter([indiv.rank for indiv in self.survivors[j]],
                            [getattr(indiv, this_attr)[i] for indiv in self.survivors[j]], c=colors[j], alpha=0.5)
            plt.title(this_attr+': '+objective_name)
        for i, feature_name in enumerate(self.feature_names):
            this_attr = 'features'
            plt.figure()
            for j, population in enumerate(self.history):
                plt.scatter([indiv.rank for indiv in population],
                            [getattr(indiv, this_attr)[i] for indiv in population],
                            c=colors[j], alpha=0.1)
                plt.scatter([indiv.rank for indiv in self.survivors[j]],
                            [getattr(indiv, this_attr)[i] for indiv in self.survivors[j]], c=colors[j], alpha=0.5)
            plt.title(this_attr+': '+feature_name)
        plt.show()
        plt.close()

    def nan2None(self, attr):
        """
        Convert from numpy nan to Python None.
        :param attr: any
        :return: any
        """
        if np.isnan(attr):
            return None
        else:
            return attr

    def None2nan(self, attr):
        """
        Convert from Python None to numpy nan.
        :param attr: any
        :return: any
        """
        if attr is None:
            return np.nan
        else:
            return attr

    def save(self, file_path):
        """

        :param file_path: str
        """
        with h5py.File(file_path, 'w') as f:
            f.attrs['param_names'] = self.param_names
            f.attrs['feature_names'] = self.feature_names
            f.attrs['objective_names'] = self.objective_names
            for pop_index in range(len(self.history)):
                f.create_group(str(pop_index))
                for group_name, population in zip(['population', 'survivors'],
                                                  [self.history[pop_index], self.survivors[pop_index]]):
                    f[str(pop_index)].create_group(group_name)
                    for ind_index, individual in enumerate(population):
                        f[str(pop_index)][group_name].create_group(str(ind_index))
                        f[str(pop_index)][group_name][str(ind_index)].attrs['energy'] = self.None2nan(individual.energy)
                        f[str(pop_index)][group_name][str(ind_index)].attrs['rank'] = self.None2nan(individual.rank)
                        f[str(pop_index)][group_name][str(ind_index)].attrs['distance'] = \
                            self.None2nan(individual.distance)
                        f[str(pop_index)][group_name][str(ind_index)].attrs['fitness'] = \
                            self.None2nan(individual.fitness)
                        f[str(pop_index)][group_name][str(ind_index)].attrs['survivor'] = \
                            self.None2nan(individual.survivor)
                        f[str(pop_index)][group_name][str(ind_index)].create_dataset('x', data=individual.x,
                                                                                     compression='gzip',
                                                                                     compression_opts=9)
                        f[str(pop_index)][group_name][str(ind_index)].create_dataset('features',
                                                                                     data=individual.features,
                                                                                     compression='gzip',
                                                                                     compression_opts=9)
                        f[str(pop_index)][group_name][str(ind_index)].create_dataset('objectives',
                                                                                     data=individual.objectives,
                                                                                     compression='gzip',
                                                                                     compression_opts=9)
        print 'PopulationStorage: saved %i generations to file: %s' % (len(self.history), file_path)

    def load(self, file_path):
        """

        :param file_path: str
        """
        with h5py.File(file_path, 'r') as f:
            self.param_names = f.attrs['param_names']
            self.feature_names = f.attrs['feature_names']
            self.objective_names = f.attrs['objective_names']
            for pop_index in range(len(f)):
                population_list, survivors_list = [], []
                for group_name, population in zip(['population', 'survivors'], [population_list, survivors_list]):
                    group = f[str(pop_index)][group_name]
                    for ind_index in range(len(group)):
                        ind_data = group[str(ind_index)]
                        individual = Individual(ind_data['x'][:])
                        individual.features = ind_data['features'][:]
                        individual.objectives = ind_data['objectives'][:]
                        individual.energy = self.nan2None(ind_data.attrs['energy'])
                        individual.rank = self.nan2None(ind_data.attrs['rank'])
                        individual.distance = self.nan2None(ind_data.attrs['distance'])
                        individual.fitness = self.nan2None(ind_data.attrs['fitness'])
                        individual.survivor = self.nan2None(ind_data.attrs['survivor'])
                        population.append(individual)
                self.history.append(population_list)
                self.survivors.append(survivors_list)
        print 'PopulationStorage: loaded %i generations to file: %s' % (len(self.history), file_path)


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


def sort_by_absolute_energy(population):
    """
    Modifies in place the energy attribute of each Individual in the population. Sorts the population by absolute
    energy (sum of all objectives). Returns the sorted population.
    :param population: list of :class:'Individual'
    """
    indexes = range(len(population))
    energy_vals = []
    for individual in population:
        this_energy = np.sum(individual.objectives)
        individual.energy = this_energy
        energy_vals.append(this_energy)
    indexes.sort(key=energy_vals.__getitem__)
    population = map(population.__getitem__, indexes)
    return population


def sort_by_relative_energy(population):
    """
    Modifies in place the energy attribute of each Individual in the population. Sorts the population by relative
    energy (sum of all normalized objectives). Returns the sorted population.
    :param population: list of :class:'Individual'
    """
    pop_size = len(population)
    num_objectives = [len(individual.objectives) for individual in population if individual.objectives is not None]
    if len(num_objectives) < pop_size:
        raise Exception('sort_by_relative_energy: objectives have not been stored for all Individuals in population')
    num_objectives = max(num_objectives)
    for individual in population:
        individual.energy = 0
    for m in range(num_objectives):
        objective_vals = [individual.objectives[m] for individual in population]
        objective_min = min(objective_vals)
        objective_max = max(objective_vals)
        if objective_min != objective_max:
            objective_vals = np.subtract(objective_vals, objective_min)
            objective_vals /= objective_max - objective_min
            for energy, individual in zip(objective_vals, population):
                individual.energy += energy
    indexes = range(pop_size)
    energy_vals = [individual.energy for individual in population]
    indexes.sort(key=energy_vals.__getitem__)
    population = map(population.__getitem__, indexes)
    return population


def assign_rank_by_fitness_and_energy(population):
    """
    Modifies in place the rank attributes of each Individual in the population. Within each group of Individuals with
    equivalent fitness, sorts by total energy (sum of all objectives).
    :param population: list of :class:'Individual'
    """
    pop_size = len(population)
    fitness_vals = [individual.fitness for individual in population]
    if len(fitness_vals) < pop_size:
        raise Exception('assign_rank_by_fitness_and_energy: fitness has not been stored for all Individuals in '
                        'population')
    max_fitness = max(fitness_vals)
    new_population = []
    for fitness in range(max_fitness + 1):
        new_front = [individual for individual in population if individual.fitness == fitness]
        new_sorted_front = sort_by_relative_energy(new_front)
        new_population.extend(new_sorted_front)
    # now that population is sorted, assign rank to Individuals
    for rank, individual in enumerate(new_population):
        individual.rank = rank


def assign_rank_by_fitness_and_crowding_distance(population):
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
        raise Exception('assign_rank_by_fitness_and_crowding_distance: fitness has not been stored for all Individuals '
                        'in population')
    max_fitness = max(fitness_vals)
    if max_fitness > 0:
        new_population = []
        for fitness in range(max_fitness + 1):
            new_front = [individual for individual in population if individual.fitness == fitness]
            new_sorted_front = sort_by_crowding_distance(new_front)
            new_population.extend(new_sorted_front)
    else:
        new_population = sort_by_relative_energy(population)
    # now that population is sorted, assign rank to Individuals
    for rank, individual in enumerate(new_population):
        individual.rank = rank


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


def evaluate_bgen(population, disp=False):
    """
    Modifies in place the fitness, energy and rank attributes of each Individual in the population.
    :param population: list of :class:'Individual'
    :param disp: bool
    """
    assign_fitness_by_dominance(population)
    assign_rank_by_fitness_and_energy(population)


def evaluate_random(population, disp):
    """
    Modifies in place the rank attribute of each Individual in the population.
    :param population: list of :class:'Individual'
    :param disp: bool
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
    parameter arrays for parallel evaluation. Features fitness-based pruning and adaptive reduction of step_size every
    iteration. Each iteration consists of path_length number of generations without pruning.
    """
    def __init__(self, x0, param_names, feature_names, objective_names, pop_size, bounds=None, take_step=None,
                 evaluate=None, seed=None, max_iter=None, max_gens=None, path_length=1, adaptive_step_factor=0.9,
                 niter_success=None, survival_rate=0.1, disp=False):
        """

        :param x0: array
        :param param_names: list of str
        :param feature_names: list of str
        :param objective_names: list of str
        :param pop_size: int
        :param bounds: list of tuple
        :param take_step: callable
        :param evaluate: callable
        :param seed: int or :class:'np.random.RandomState'
        :param max_iter: int
        :param max_gens: int
        :param path_length: int
        :param adaptive_step_factor: float in [0., 1.]
        :param niter_success: int
        :param survival_rate: float in [0., 1.]
        :param disp: bool
        """
        self.x0 = np.array(x0)
        self.storage = PopulationStorage(param_names=param_names, feature_names=feature_names,
                                         objective_names=objective_names)
        self.pop_size = pop_size
        if evaluate is None:
            self._evaluate = evaluate_bgen
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
        self.path_length = path_length
        if max_iter is None:
            if max_gens is None:
                self.max_gens = self.path_length * 30
            else:
                self.max_gens = max_gens
        else:
            if max_gens is None:
                self.max_gens = self.path_length * max_iter
            else:
                self.max_gens = min(max_gens, self.path_length * max_iter)
        self.adaptive_step_factor = adaptive_step_factor
        if niter_success is None:
            self.ngen_success = self.max_gens
        else:
            self.ngen_success = min(self.max_gens, self.path_length * niter_success)
        self.num_survivors = max(1, int(pop_size * survival_rate))
        self.disp = disp
        self.objectives_stored = False
        self.evaluated = False
        self.population = []
        self.survivors = []
        self.final_survivors = None

    def __call__(self):
        """
        A generator that yields a list of size pop_size of parameter arrays.
        :yields: list of :class:'Individual'
        """
        self.num_gen = 0
        self.start_time = time.time()
        self.local_time = self.start_time
        while self.num_gen < self.max_gens:
            if self.num_gen == 0:
                self.init_population()
                self.evaluated = False
            else:
                if not self.objectives_stored:
                    raise Exception('BGen: Gen %i, objectives have not been stored for all Individuals in '
                                    'population' % self.num_gen)
                if self.disp:
                    print 'BGen: Gen %i, computing features took %.2f s' % \
                          (self.num_gen - 1, time.time() - self.local_time)
                if self.num_gen % self.path_length == 0:
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
            self.local_time = time.time()
            self.num_gen += 1
            yield [individual.x for individual in self.population]
        # evaluate the final, potentially incomplete interval of generations
        if self.objectives_stored and not self.evaluated:
            if self.disp:
                print 'BGen: Gen %i, computing features took %.2f s' % (self.num_gen-1, time.time()-self.local_time)
            generations = self.num_gen % self.path_length
            if generations == 0:
                generations = self.path_length
            self.final_survivors = self.storage.get_best(self.num_survivors,
                                  generations=generations, evaluate=self._evaluate,
                                  modify=True)
            for individual in self.final_survivors:
                individual.survivor = True
            self.evaluated = True
            if self.disp:
                print 'BGen: %i generations took %.2f s' % (self.max_gens, time.time()-self.start_time)

    def update_population(self, features, objectives):
        """
        Expects a list of objective arrays to be in the same order as the list of parameter arrays yielded from the
        current generation.
        :param features: list of dict
        :param objectives: list of dict
        """
        for i, objective_dict in enumerate(objectives):
            if type(objective_dict) != dict:
                raise TypeError('BGen.update_population: objectives must be a list of dict')
            this_objectives = np.array([objective_dict[key] for key in self.storage.objective_names])
            self.population[i].objectives = this_objectives
        for i, feature_dict in enumerate(features):
            if type(feature_dict) != dict:
                raise TypeError('BGen.update_population: features must be a list of dict')
            this_features = np.array([feature_dict[key] for key in self.storage.feature_names])
            self.population[i].features = this_features
        if self.disp:
            print 'BGen: Gen %i, storing features and objectives for population size %i' % \
                  (self.num_gen-1, self.pop_size)
        self.storage.append(self.population, self.survivors)
        self.objectives_stored = True

    def evaluate(self):
        """
        Assign fitness
        """
        self._evaluate(self.population + self.survivors, disp=self.disp)
        self.evaluated = True

    def init_population(self):
        """

        """
        self.population = [Individual(self.take_step(self.x0)) for i in range(self.pop_size)]

    def step_survivors(self):
        """
        Consider the highest ranked Individuals of the previous interval of generations to be the survivors. Seed the
        next generation with steps taken from the surviving set of parameters.
        """
        self.survivors = self.storage.get_best(self.num_survivors,
                                               generations=self.path_length, evaluate=self._evaluate,
                                               modify=True)
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

    def report_best(self, generations=None):
        """
        Format and print the contents of self.final_survivors.
        """
        pass


class EGen(object):
    """
    This class is inspired by emoo (Bahl A, Stemmler MB, Herz AVM, Roth A. (2012). J Neurosci Methods)
    """

    def __init__(self, x0, param_names, feature_names, objective_names, pop_size, bounds=None, interval=None, custom=None):
        """
        self, x0, param_names, feature_names, objective_names, pop_size, bounds=None, take_step=None,
                 evaluate=None, seed=None, max_iter=None, max_gens=None, path_length=1, adaptive_step_factor=0.9,
                 niter_success=None, survival_rate=0.1, disp=False):

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

            for i in np.arange(0, len(self.population) - 1, 2):
                individual1, individual2 = population_permutation[i], population_permutation[i + 1]
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
        if (self.mpi == False):

            for individual in self.population:

                # only evaluate those that are really new!
                if individual[self.fitnesspos] == -1:

                    parameters = individual[:self.para]

                    objectives_error = self.evaluate_individual(parameters)

                    if (objectives_error != None):
                        new_population.append(np.r_[parameters, objectives_error, self.no_properties])
                else:
                    new_population.append(individual)
        else:
            # distribute the individuals among the slaves
            i = 0
            for individual in self.population:
                if individual[self.fitnesspos] == -1:
                    parameters = individual[:self.para]

                    dest = i % (self.comm.size - 1) + 1
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

        while (True):
            parameters = self.comm.recv(source=0)  # wait....

            # Does the master want the slave to shutdown?
            if (parameters == None):
                # Slave finishing...
                break

            objectives_error = self.evaluate_individual(parameters)

            # objectives_error = self.get_objectives_error(self.unnormit(parameters))
            if (objectives_error == None):
                self.comm.send(None, dest=0)
            else:
                self.comm.send([parameters, objectives_error], dest=0)

    def assign_fitness(self):
        """
        are we in a multiobjective regime, then the selection of the best individual is not trival
        and must be based on dominance, thus we determine all non dominated fronts and only use the best
        to transfer into the new generation
        """
        if (self.obj > 1):
            self.assign_rank()

            new_population = np.array([])

            maxrank = self.population[:, self.rankpos].max()

            for rank in range(0, int(maxrank) + 1):

                new_front = self.population[np.where(self.population[:, self.rankpos] == rank)]

                new_sorted_front = self.crowding_distance_sort(new_front)

                if (len(new_population) == 0):
                    new_population = new_sorted_front
                else:
                    new_population = np.r_[new_population, new_sorted_front]

            self.population = new_population

        else:
            # simple sort the objective value
            ind = np.argsort(self.population[:, self.objpos])
            self.population = self.population[ind]

        # now set the fitness, indiviauls are sorted, thus fitnes is easy to set
        fitness = range(0, len(self.population[:, 0]))
        self.population[:, -1] = fitness

    def new_generation(self):
        # the worst are at the end, let them die, if there are too many
        if (len(self.population) > self.size):
            self.population = self.population[:self.size]

    def dominates(self, p, q):

        objectives_error1 = self.population[p][self.objpos:self.objpos + self.obj]
        objectives_error2 = self.population[q][self.objpos:self.objpos + self.obj]

        diff12 = objectives_error1 - objectives_error2

        # is individdum equal or better then individdum two?
        # and at least in one objective better
        # then it dominates individuum2
        # if not it does not dominate two (which does not mean that 2 may not dominate 1)
        return (((diff12 <= 0).all()) and ((diff12 < 0).any()))

    def assign_rank(self):

        F = dict()

        P = self.population

        S = dict()
        n = dict()
        F[0] = []

        # determine how many solutions are dominated or dominate
        for p in range(len(P)):

            S[p] = []  # this is the list of solutions dominated by p
            n[p] = 0  # how many solutions are dominating p

            for q in range(len(P)):

                if self.dominates(p, q):
                    S[p].append(q)  # add q to the list of solutions dominated by p
                elif self.dominates(q, p):
                    n[p] += 1  # q dominates p, thus increase number of solutions that dominate p

            if n[p] == 0:  # no other solution dominates p

                # this is the rank column
                P[p][self.rankpos] = 0

                F[0].append(p)  # add p to the list of the first front

        # find the other non dominated fronts
        i = 0
        while len(F[i]) > 0:
            Q = []  # this will be the next front

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
            F[i] = Q  # this is the next front

    def crowding_distance_sort(self, front):

        sorted_front = front.copy()

        l = len(sorted_front[:, 0])

        sorted_front[:, self.distpos] = np.zeros_like(sorted_front[:, 0])

        for m in range(self.obj):
            ind = np.argsort(sorted_front[:, self.objpos + m])
            sorted_front = sorted_front[ind]

            # definitely keep the borders
            sorted_front[0, self.distpos] += 1000000000000000.
            sorted_front[-1, self.distpos] += 1000000000000000.

            fm_min = sorted_front[0, self.objpos + m]
            fm_max = sorted_front[-1, self.objpos + m]

            if fm_min != fm_max:
                for i in range(1, l - 1):
                    sorted_front[i, self.distpos] += (sorted_front[i + 1, self.objpos + m] - sorted_front[
                        i - 1, self.objpos + m]) / (fm_max - fm_min)

        ind = np.argsort(sorted_front[:, self.distpos])
        sorted_front = sorted_front[ind]
        sorted_front = sorted_front[-1 - np.arange(len(sorted_front))]

        return sorted_front
