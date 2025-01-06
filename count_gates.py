import os.path
import numpy as np
import VQE
from codon_optimization import CodonOptimization
from one_hot_codon_optimization import CodonOptimizationOneHot
from tunable_parameters import TunableParameters
import util, csv, time

def write_data(name, data):
    with open(name, 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(data)


def count_gates(protein_sequence_list, encode_type):
    results = [["amino_sequence", "qubit_len", "Depth", "Gate_counts"]]
    index = 1
    count = len(protein_sequence_list)
    for sequence in protein_sequence_list:
        tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
        if encode_type == "one_hot":
            codonOpt = CodonOptimizationOneHot(sequence, tune_parameters)

        else:
            codonOpt = CodonOptimization(sequence, tune_parameters)
        print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
        qubit_op = codonOpt.create_qubit_op()
        ansatz = VQE.get_ansatz(qubit_op).decompose()
        gate_count = ansatz.count_ops()
        total_gate = sum(gate_count.values())
        index += 1
        print("Gate counts:", total_gate)
        results.append([sequence, codonOpt.qubit_len, ansatz.depth(), total_gate])
    return results


def spilt_sequence(protein_sequence, split_len):
    sequence_list = []

    index = 0
    while index < len(protein_sequence):
        sequence_list.append(protein_sequence[index:index + split_len])
        index += split_len

    return sequence_list


def write_gates(protein_sequence, split_range, encode_type):
    # statistic the number of gates with dense encoding for fragments with a length from 3 to 19
    for split_len in split_range:
        sequence_list = spilt_sequence(protein_sequence, split_len)
        output_data = count_gates(sequence_list, encode_type)
        # name = time.strftime("%y-%m-%d_%H-%M-%S")
        path = f"output/{split_len}/Gate"
        if not os.path.exists(path):
            os.makedirs(path)
        write_data(f"{path}/{encode_type}_P0DTC2_Gate.csv", output_data)


def get_mean_max(file_name, target):
    value = []
    with open(file_name, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Sum the values in the target column (convert to float for numerical data)
            value.append(float(row[target]))

    return np.mean(value), max(value)


def statistic_gates(encode_type, split_range):
    result = [["fragment_len", "Mean_Depth", "Max_Depth", "Mean_Gate", "Max_Depth"]]
    for split_len in split_range:
        file_name = f"output/{split_len}/Gate/{encode_type}_P0DTC2_Gate.csv"
        if not os.path.exists(file_name):
            print(f"{file_name} is not exist!")
        else:
            mean_gate, max_gate = get_mean_max(file_name, "Gate_counts")
            mean_depth, max_depth = get_mean_max(file_name, "Depth")
            result.append([split_len, np.ceil(mean_depth), max_depth, np.ceil(mean_gate), max_gate])
    write_data(f"output/{encode_type}_depth_gate.csv", result)


if __name__ == '__main__':
    protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
    # This is to statistic gates and qubits for the paper A resource-efficient variational quantum algorithm for mRNA
    # codon optimization
    split_range = range(6, 20)
    encode_type = "one_hot"
    # write_gates(protein_sequence, split_range, encode_type)
    statistic_gates(encode_type, split_range)
    # statistic_qubits(protein_sequence, split_range)
    encode_type = "dense"
    # write_gates(protein_sequence, split_range, encode_type)
    statistic_gates(encode_type, split_range)