import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import util
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



class GA:

    def __init__(self, protein_sequence, tune_parameters):
        self.protein_sequence = protein_sequence
        self.tune_parameters = tune_parameters
        self.codonOpt = CodonOptimization(protein_sequence, tune_parameters)
        self.qubit_op = self.codonOpt.create_qubit_op()
        self.qubit_len = self.codonOpt.qubit_len

    def exact(self, X):
        opt = 0

        for op in self.qubit_op.primitive:
            item = 1
            loc = find_all(op.paulis[0])
            for index in loc:
                item *= X[index]
            opt += item * op.coeffs[0].real
        return opt

    def compute(self):
        energy_min = 10**5
        sequence = "0" * self.qubit_len

        for i in range(2**self.qubit_len):
            binary = util.int_to_binary(i, self.qubit_len)
            z = [1 if i == "0" else -1 for i in binary]
            energy = self.exact(z)
            if energy_min < energy:
                energy_min = energy
                sequence = binary
            print(binary, energy)
        print(f"The optimal sequence {sequence}, value {energy_min}")

    def execute(self):
        varbound = np.array([[0, 1]] * self.qubit_len)
        model = ga(function=self.exact, dimension=self.qubit_len, variable_type='bool', variable_boundaries=varbound)
        model.run()

tune_parameters = TunableParameters(0.3, 0.15*6, 1, 13000)
protein_sequence = "GSK"

ga = GA(protein_sequence, tune_parameters)
ga.execute()
# protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
# split_len = 8
# sequence_list = util.spilt_sequence(protein_sequence, split_len)
# for sequence in sequence_list:
#     ga = GA(sequence, tune_parameters)
#     ga.execute()
#     break