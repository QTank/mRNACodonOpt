import VQE
from codon_optimization import CodonOptimization
from tunable_parameters import TunableParameters
import util
import csv
import numpy as np
import time


def write_data(name, data):
    with open(name, 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["amino_sequence", "dense_encoding", "qubit_len", "mRNA_sequence", "value", "time", "Gate_counts"])

        writer.writerows(data)

def optimization(protein_sequence):
    results = []
    index = 1
    count = len(protein_sequence)
    results.append(["amino_sequence", "dense_encoding", "qubit_len", "mRNA_sequence", "value", "time"])
    for sequence in protein_sequence:
        start_time = time.time()
        tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
        codonOpt = CodonOptimization(sequence, tune_parameters)
        print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
        qubit_op = codonOpt.create_qubit_op()
        count_pauli_z = util.count_pauli_z(qubit_op)
        bitstring, value = VQE.get_min(qubit_op)
        # bitstring, value = VQE.test(qubit_op)
        mRNA_sequence = util.decode_bitstring(sequence, bitstring)

        executing_time = round(time.time() - start_time, 2)
        print(f"The bit string is {bitstring}, value is {value.real}")
        print(f"The execute time on {sequence} is {executing_time} s")
        index += 1
        results.append([sequence, bitstring, len(bitstring), mRNA_sequence, round(value.real, 3), executing_time, count_pauli_z])

    return results

def count_gates(protein_sequence):
    results = []
    index = 1
    count = len(protein_sequence)
    results.append(["amino_sequence", "qubit_len", "Gate_counts"])
    for sequence in protein_sequence:
        tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)

        codonOpt = CodonOptimization(sequence, tune_parameters)
        print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
        qubit_op = codonOpt.create_qubit_op()
        count_pauli_z = util.count_pauli_z(qubit_op)
        index += 1
        results.append([sequence, codonOpt.qubit_len, count_pauli_z])
    return results

def spilt_sequence(protein_sequence, split_len):
    sequence_list = []

    index = 0
    while index < len(protein_sequence):
        sequence_list.append(protein_sequence[index:index+split_len])
        index += split_len

    return sequence_list

protein_sequence = "KIILFLALITLATCELYHYQECVRGTTVLLKEPCSSGTYEGNSPFHPLADNKFALTCFSTQFAFACPDGVKHVYQLRARSVSPKLFIRQEEVQELYSPIFLIVAAIVFITLCFTLKRKTE"
protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

for i in range(19, 20):
    split_len = i
    if split_len == 8:
        continue

    sequence_list = spilt_sequence(protein_sequence, split_len)
    name = time.strftime("%m-%d_%H-%M-%S")
# output_data = optimization(sequence_list)
    output_data = count_gates(sequence_list)
# write_data(f"output/{split_len}/P0DTC2_{name}.csv", output_data)
    write_data(f"output/{split_len}/Gate/dense_P0DTC2_{name}_Gate.csv", output_data)

# sequence = "HAIHVSGT"
# parameters = TunableParameters(0.3, 0.15 * 6, 1, 13000)
# codonOpt = CodonOptimization(sequence, parameters)
# qubit_op = codonOpt.create_qubit_op()
# bitstring, value = VQE.get_min(qubit_op)
# mRNA_sequence = util.decode_bitstring(sequence, bitstring)
# print(bitstring, value)
# print(mRNA_sequence)
