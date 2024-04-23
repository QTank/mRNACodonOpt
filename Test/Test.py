import python_codon_tables as pct
from tunable_parameters import TunableParameters
import numpy as np

import util, time
from codon_optimization import CodonOptimization

class Test:
    def __init__(self, protein, binary_str, tunable_parameters, table_name='e_coli_316407'):
        self.protein_sequence = protein
        self.binary_str = binary_str
        self.tunable_parameters = tunable_parameters
        self.codon_table = pct.get_codons_table(table_name)
        self.constant_for_usage = 0.1
        self.target_gc = 0.5
        self.nucleotides_len = 3 * len(protein)
        self.amino_codon_list = util.decode_bitstring_list(protein, binary_str) if self.judge_binary_str(binary_str) else []

    def get_energy(self):
        if self.judge_binary_str(self.binary_str):
            return self.get_h_usage() + self.get_target_gc() + self.get_repeated()
        return 130000

    def decode_codon(self, qbit):
        start_index = 0

        amino_codon_list = []
        for amino in self.protein_sequence:
            codon_list = list(self.codon_table[amino].keys())
            codon_encoding_len = round(np.log2(len(codon_list)))
            if codon_encoding_len == 0:
                amino_codon_list.append(codon_list[0])
            else:
                loc = int(qbit[start_index:start_index+codon_encoding_len], 2)
                start_index += codon_encoding_len
                amino_codon_list.append(codon_list[loc])

        return amino_codon_list

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

        try:
            for amino in self.protein_sequence:
                codon_list = list(self.codon_table[amino].keys())
                codon_len = round(np.log2(len(codon_list)))
                if codon_len == 0:
                    continue
                loc = int(qbit[start_index:start_index+codon_len], 2)
                start_index += codon_len
                if loc >= len(codon_list):
                    # print(f"The binary string {qbit} is invalid!")
                    return False
            return True
        except:
            return False


def iterate_all_value(amino_sequence, encoding_len):
    opt_value = 1000000
    bitstring = ""
    for i in range(2 ** encoding_len - 1):
        binary_str = util.int_to_binary(i, encoding_len)
        if binary_str != "00000010100":
            continue
        test = Test(amino_sequence, binary_str, tunable_parameters)

        if test.judge_binary_str(binary_str):
            v = test.get_energy()
            # print(f"the bitsring {binary_str}, value {v}, {util.decode_bitstring(amino_sequence, binary_str)}")
            if v < opt_value:
                opt_value, bitstring = v, binary_str
    print(f"The opt string: {bitstring}, value: {opt_value}")
    return bitstring, opt_value


tunable_parameters = TunableParameters(0.3, 0.15 * 6, 1, 13000)
split_len = 8
name = time.strftime("%m-%d_%H-%M-%S")
value_type = "exact"
name = f"output/{split_len}/{value_type}/P0DTC2_{name}.csv"

protein_sequence_covid = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
protein_sequence_list = util.spilt_sequence(protein_sequence_covid, split_len)

result = []
for sequence in protein_sequence_list[10:]:

    sequence = "FFSNVTWF"
    start_time = time.time()
    codon = CodonOptimization(sequence, tunable_parameters)
    encoding_len = codon.qubit_len

    bitstring, opt_value = iterate_all_value(sequence, encoding_len)

    mrna_sequence = util.decode_bitstring(sequence, bitstring)
    running_time = time.time() - start_time
    result.append([sequence, bitstring, encoding_len, mrna_sequence, opt_value, running_time])
    # util.write_data(name, result)
    break