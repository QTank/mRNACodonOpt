import os.path
import VQE
from codon_optimization import CodonOptimization
from tunable_parameters import TunableParameters
import util, csv, time
import numpy as np
from one_hot_codon_optimization import CodonOptimizationOneHot


def write_data(name, data):
    with open(name, 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(
            ["amino_sequence", "dense_encoding", "qubit_len", "mRNA_sequence", "value", "time", "Gate_counts"])

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
        mRNA_sequence = util.decode_bitstring(sequence, bitstring)

        executing_time = round(time.time() - start_time, 2)
        print(f"The bit string is {bitstring}, value is {value.real}")
        print(f"The execute time on {sequence} is {executing_time} s")
        index += 1
        results.append(
            [sequence, bitstring, len(bitstring), mRNA_sequence, round(value.real, 3), executing_time, count_pauli_z])

    return results


def count_gates(protein_sequence_list):
    results = []
    index = 1
    count = len(protein_sequence_list)
    results.append(["amino_sequence", "qubit_len", "Gate_counts"])
    for sequence in protein_sequence_list:
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
        sequence_list.append(protein_sequence[index:index + split_len])
        index += split_len

    return sequence_list


def statistic_gates(protein_sequence, split_range):
    # statistic the number of gates with dense encoding for fragments with a length from 3 to 19
    for split_len in split_range:
        sequence_list = spilt_sequence(protein_sequence, split_len)
        output_data = count_gates(sequence_list)
        name = time.strftime("%m-%d_%H-%M-%S")
        path = f"../output/{split_len}/Gate"
        if not os.path.exists(path):
            os.makedirs(path)
        write_data(f"{path}/dense_P0DTC2_{name}_Gate.csv", output_data)


def execute_optimization(protein_sequence, split_range):
    # execute the quantum method for mRNA codon optimization for fragments with a length from 3 to 8.
    # It takes more time when the length of fragments is larger than 8
    for split_len in split_range:
        sequence_list = spilt_sequence(protein_sequence, split_len)
        output_data = optimization(sequence_list)
        name = time.strftime("%m-%d_%H-%M-%S")
        path = f"../output/{split_len}/dense"
        if not os.path.exists(path):
            os.makedirs(path)
        write_data(f"/P0DTC2_{name}.csv", output_data)


def statistic_qubits(protein_sequence, split_range):
    # statistic the required number of qubits for each fragment,
    # the maximum qubits of each fragment as result,
    # the mean qubits of each fragment as result,
    max_qubits = [["fragment_len", "dense", "one_hot"]]
    mean_qubits = [["fragment_len", "dense", "one_hot"]]
    for split_len in split_range:
        sequence_list = spilt_sequence(protein_sequence, split_len)
        qubits_dense_list = util.count_qubits(sequence_list, "dense")
        qubits_one_hot_list = util.count_qubits(sequence_list, "one_hot")
        max_qubits.append([split_len, max(qubits_dense_list), max(qubits_one_hot_list)])
        mean_qubits.append([split_len, round(np.mean(qubits_dense_list)), round(np.mean(qubits_one_hot_list))])
    util.write_data(f"../output/statistic_qubits.csv", max_qubits)
    util.write_data(f"../output/statistic_mean_qubits.csv", mean_qubits)


def solution(sequence):
    # test whether the quantum method can work for the protein sequence "HAIHVSGT"
    parameters = TunableParameters(0.3, 0.15 * 6, 1, 13000)
    codonOpt = CodonOptimization(sequence, parameters)
    qubit_op = codonOpt.create_qubit_op()
    bitstring, value = VQE.get_min(qubit_op)
    mRNA_sequence = util.decode_bitstring(sequence, bitstring)
    print(bitstring, value)
    print(mRNA_sequence)


if __name__ == '__main__':
    protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
    split_range = range(6, 20)
    # statistic_gates(protein_sequence, split_range)
    # statistic_qubits(protein_sequence, split_range)
    # execute_optimization(protein_sequence, range(3, 9))
    sequence = "HAIHVSGT"
    solution(sequence)

