import python_codon_tables as pct
from tunable_parameters import TunableParameters
import numpy as np

import util
from codon_optimization import CodonOptimization

class Test:
    def __init__(self, protein, binary_str, tunable_parameters, table_name='e_coli_316407'):
        self.protein_sequence = protein
        self.binary_str = binary_str
        self.tunable_parameters = tunable_parameters
        self.codon_table = pct.get_codons_table(table_name)
        self.constant_for_usage = 0.1
        self.target_gc = 0.5
        self.nucleotides_len = 3 * len(protein_sequence)
        self.amino_codon_list = self.decode_codon(binary_str) if self.judge_binary_str(binary_str) else {}

    def get_energy(self):
        return self.get_h_usage() + self.get_target_gc() + self.get_repeated()

    def decode_codon(self, qbit):
        start_index = 0

        mRNA_codon_list = []
        for amino in self.protein_sequence:
            codon_list = list(self.codon_table[amino].keys())
            encoding_len = round(np.log2(len(codon_list)))
            loc = int(qbit[start_index:start_index+encoding_len], 2)
            start_index += encoding_len

            codon_sequence = codon_list[loc]
            mRNA_codon_list.append(codon_sequence)

        return mRNA_codon_list

    def get_h_usage(self):
        h_usage = 0

        for index in range(len(self.protein_sequence)):
            amino = self.protein_sequence[index]
            codon_sequence = self.amino_codon_list[index]
            h_usage += round(-np.log10(self.codon_table[amino][codon_sequence] + self.constant_for_usage), 3)
        return h_usage * self.tunable_parameters.get_usage() * (self.nucleotides_len ** 2)

    def get_target_gc(self):
        h_gc = 0

        for codon_sequence in self.amino_codon_list:
            h_gc += util.get_gc_count(codon_sequence)

        return ((h_gc - self.target_gc * self.nucleotides_len) ** 2) * self.tunable_parameters.get_target_gc()

    def get_repeated(self):
        h_repeated = 0

        for curr_index in range(1, len(self.protein_sequence)):
            pre_codon = self.amino_codon_list[curr_index-1]
            curr_codon = self.amino_codon_list[curr_index]
            h_repeated += util.get_neighbour_repetition(pre_codon, curr_codon)

        return h_repeated * self.tunable_parameters.get_repeated_nucleotides() * (self.nucleotides_len ** 2)

    def judge_binary_str(self, qbit):
        start_index = 0

        for amino in self.protein_sequence:
            codon_list = list(self.codon_table[amino].keys())
            encoding_len = round(np.log2(len(codon_list)))
            loc = int(qbit[start_index:start_index+encoding_len], 2)
            start_index += encoding_len
            if loc >= len(codon_list):
                print(f"The binary string {qbit} is invalid!")
                return False
        return True


def get_opt_bitstring(protein_sequence):
    tunable_parameters = TunableParameters(0.3, 0.15 * 6, 1, 13000)
    protein_sequence = "FVF"
    codon = CodonOptimization(protein_sequence, tunable_parameters)
    bitstring, opt_value = "", 13000
    encoding_len = codon.qubit_len

    for i in range(2**encoding_len - 1):
        binary_str = util.int_to_binary(i, encoding_len)

        test = Test(protein_sequence, binary_str, tunable_parameters)

        if test.judge_binary_str(binary_str):
            v = test.get_energy()
            # print(f"the bitsring {binary_str}, value {v}, decoding{util.decode_bitstring(protein_sequence, binary_str)}")
            if v < opt_value:
                opt_value, bitstring = v, binary_str
    print(f"The opt string: {bitstring}, value: {opt_value}")
    return binary_str, opt_value


def write_opt_bitstring():
    start_time = time.time()
    tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
    codonOpt = CodonOptimizationOneHot(sequence, tune_parameters)
    print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
    qubit_op = codonOpt.create_qubit_op()
    bitstring, value = VQE.get_min(qubit_op)
    mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring)

    executing_time = round(time.time() - start_time, 2)
    print(f"The bit string is {bitstring}, value is {value.real}")
    print(f"The execute time on {sequence} is {executing_time} s")

