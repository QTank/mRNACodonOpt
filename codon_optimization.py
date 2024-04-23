import numpy as np
import python_codon_tables as pct
from codon import Codon
from tunable_parameters import TunableParameters
import util

class CodonOptimization:

    def __init__(self, protein_sequence: str, parameters: TunableParameters, table_name="e_coli_316407"):
        self.protein_sequence = protein_sequence
        self.nucleotides_len = 3 * len(protein_sequence)
        self.tunable_parameters = parameters

        self.codon_table = pct.get_codons_table(table_name)
        self.qubit_len = self.calculate_qubit_len()
        self.codon_list = self.create_codon_list()

        self.constant_for_usage = 0.1
        self.target_gc = 0.5

    def calculate_qubit_len(self):
        count = 0
        for amino in self.protein_sequence:
            count += round(np.log2(len(self.codon_table[amino])))

        return count

    def create_codon_list(self):

        codon_list = []
        loc_index = 0
        for index in range(len(self.protein_sequence)):
            codon = Codon(self.protein_sequence[index], loc_index, self.qubit_len)
            loc_index += codon.encoding_qubit_len
            codon_list.append(codon)
        return codon_list

    def create_usage(self):
        h_usage = 0

        for codon in self.codon_list:
            if codon.amino in ["W", "M"]:
                continue

            for sequence in codon.indicator_dict:
                indicator_frequency = codon.indicator_dict[sequence]
                h_usage += indicator_frequency[0].reduce() * round(-np.log10(indicator_frequency[1] + self.constant_for_usage), 3)

        tune = self.tunable_parameters.get_usage() * (self.nucleotides_len ** 2)
        return h_usage * tune

    def create_target_gc(self):
        h_gc = 0

        for codon in self.codon_list:
            for sequence in codon.indicator_dict:
                indicator_frequency = codon.indicator_dict[sequence]
                h_gc += indicator_frequency[0] * util.get_gc_count(sequence)

        return ((h_gc - util.build_full_identity(self.qubit_len) * self.target_gc * self.nucleotides_len) ** 2) * self.tunable_parameters.get_target_gc()

    def create_repeated_nucleotides(self):
        h_repeated = 0

        for current_index in range(1, len(self.protein_sequence)):
            pre_index = current_index - 1

            for pre_codon_sequence in self.codon_list[pre_index].indicator_dict:
                pre_indicator_frequency = self.codon_list[pre_index].indicator_dict[pre_codon_sequence]
                for cur_codon_sequence in self.codon_list[current_index].indicator_dict:
                    cur_indicator_frequency = self.codon_list[current_index].indicator_dict[cur_codon_sequence]
                    indicator_qubit = pre_indicator_frequency[0] @ cur_indicator_frequency[0]
                    h_repeated += indicator_qubit * util.get_neighbour_repetition(pre_codon_sequence, cur_codon_sequence)

        tune = self.tunable_parameters.get_repeated_nucleotides() * (self.nucleotides_len ** 2)
        return h_repeated * tune

    def create_redundant_encoding(self):
        redundant_list = []
        for codon in self.codon_list:
            redundant_list.extend(codon.redundant_indicator_list)

        return self.tunable_parameters.get_redundant_encoding() * sum(redundant_list)

    def create_qubit_op(self):
        h = self.create_usage() + self.create_target_gc() + self.create_repeated_nucleotides() + self.create_redundant_encoding()

        # h = self.create_redundant_encoding()
        return h.reduce()

