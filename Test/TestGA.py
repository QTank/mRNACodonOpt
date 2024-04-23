import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import util, time
from tunable_parameters import TunableParameters
from codon_optimization import CodonOptimization

def find_all(search):
    loc = []
    search = str(search)
    for i in range(len(search)):
        if search[i] == "Z":
            loc.append(i)

    return loc


def f(X):
    Z = [2 * (0.5 + (-1)**X[i]) for i in range(len(X))]
    opt = Z[0] - 0.2 * Z[1] + 0.3 * Z[2] - 0.7 * Z[0] * Z[1]
    return opt





def exact(X):
    opt = 0
    for op in qubit_op.primitive:
        item = 1
        loc = find_all(op.paulis[0])
        for index in loc:
            item *= X[index]
        opt += item * op.coeffs[0].real
    return opt

def compute():
    energy = []
    for i in range(2**qubit_len):
        binary = util.int_to_binary(i, qubit_len)

        z = [1 if i == "0" else -1 for i in binary]
        energy.append(exact(z))

    min = energy[0]
    sequence = "0" * qubit_len
    for i in range(1, len(energy)):
        if energy[i] < min:
            min = energy[i]
            sequence = util.int_to_binary(i, qubit_len)
        print(util.int_to_binary(i, qubit_len), energy[i])
    print(min, sequence)


# compute()
protein_sequence = "GSK"
for protein_sequence in ["GSK", "GFV"]:
    start_time = time.time()
    tune_parameters = TunableParameters(0.3, 0.15*6, 1, 13000)
    codonOpt = CodonOptimization(protein_sequence, tune_parameters)
    qubit_op = codonOpt.create_qubit_op()
    qubit_len = codonOpt.qubit_len
    qubit_len = codonOpt.qubit_len
    varbound = np.array([[0,1]]*qubit_len)
    algorithm_param = {'max_num_iteration': None, 'population_size': 100,  # None
                       'mutation_probability': 0.1, 'elit_ratio': 0.1,
                       'crossover_probability': 0.5, 'parents_portion': 0.2,
                       'crossover_type': 'uniform', 'max_iteration_without_improv': 100}


    model=ga(function=exact,dimension=qubit_len,variable_type='bool',variable_boundaries=varbound, algorithm_parameters=algorithm_param)
    model.run()

    executing_time = time.time() - start_time
    print(f"The executing time is {executing_time}")
    print(model.output_dict)