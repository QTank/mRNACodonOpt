import pygad
import numpy
import itertools


def mock_function(ga_instance, solution, sol_idx):
    return -1 * (solution[0] ** 2 + solution[1] ** 2)

range = [.1, 1.]
param_values = numpy.linspace(range[0], range[1], 10).tolist()
initial_population = list(itertools.product(param_values, param_values))

## Remove duplicates from the initial population.
for idx, sol in enumerate(initial_population):
    sol = list(sol)
    if sol[0] == sol[1]:
        if sol[1] < 0.101:
            sol[1] = sol[1] + 0.001
        else:
            sol[1] = sol[1] - 0.001
    initial_population[idx] = sol.copy()

ga = pygad.GA(num_generations=100,
              num_parents_mating=10,
              gene_type=float,
              gene_space={'low': range[0], 'high': range[1]},
              fitness_func=mock_function,
              initial_population=initial_population,
              mutation_probability=0.2,
              allow_duplicate_genes=False,
              save_solutions=True,
              suppress_warnings=True,
              mutation_by_replacement=True
              )
ga.run()

new_points = ga.population
unique_points = set([tuple(x) for x in new_points])
print(len(initial_population), len(new_points), len(unique_points))
print(f"Best solution is {ga.best_solution()[0]}.")
for point in ga.solutions:
    for x in point:
        if not range[0] <= x <= range[1]:
            print(f"Point {point} is outside gene_space.")

print("Maximum gene value found", numpy.max(ga.solutions))
print("Maximum gene value found", numpy.min(ga.solutions))