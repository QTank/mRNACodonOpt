import VQELib as VQE
from codon_optimization import CodonOptimization
from tunable_parameters import TunableParameters
import util
import random


def solution(sequence):
    # test whether the quantum method can work for the protein sequence "HAIHVSGT"
    parameters = TunableParameters(0.3, 0.15 * 6, 1, 13000)
    codon_opt = CodonOptimization(sequence, parameters)
    qubit_op = codon_opt.create_qubit_op()

    print("Starting the VQE (Variational Quantum Eigensolver) on a Quantum Simulator...")
    print("Waiting...\n")
    bitstring, minimal_value = VQE.get_min(qubit_op)
    final_mRNA_sequence = util.decode_bitstring(sequence, bitstring)

    print("Results:")
    print(f"  - Input Amino Acid Sequence: {amino_acid_sequence}")
    print(f"  - The required number of qubit: {codon_opt.qubit_len}")
    print(f"  - Minimal Value Achieved: {minimal_value.real:.2f}")
    print(f"  - Final Binary String Representation: {bitstring}")
    print(f"  - Optimized mRNA Codon Sequence: {final_mRNA_sequence}\n")
    print("Optimized Codon Mapping:")
    print("  Each amino acid is mapped to its optimized codon:")
    final_string_index = 0
    for index in range(len(amino_acid_sequence)):
        amino_acid = amino_acid_sequence[index]
        print(f"  {amino_acid:} → {bitstring[final_string_index:final_string_index + codon_opt.codon_list[index].encoding_qubit_len]} → {final_mRNA_sequence[index * 3: (index+1) * 3]}")
        final_string_index += codon_opt.codon_list[index].encoding_qubit_len
    return final_mRNA_sequence


def select_case(data):
    # Open the file and read all lines
    with open("../data/dataset", 'r') as file:
        lines = file.readlines()

    # Randomly choose one line
    return random.choice(lines).strip()


if __name__ == '__main__':
    #Test sequence HAIHVSGT
    amino_acid_sequence = "HAIHVSGT"
    # amino_acid_sequence = select_case("dataset")
    print(f"Quantum Optimization of mRNA Codons using VQE for the amino acid sequence {amino_acid_sequence}")
    final_sequence = solution(amino_acid_sequence)



