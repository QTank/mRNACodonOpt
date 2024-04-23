from one_hot_codon import CodonOneHot
import time, csv
from tunable_parameters import TunableParameters

from one_hot_codon_optimization import CodonOptimizationOneHot
import VQE
import util


def write_data(name, protein_sequence):
    with open(name, 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["amino_sequence", "dense_encoding", "qubit_len", "mRNA_sequence", "value", "time"])

        index = 1
        count = len(protein_sequence)
        for sequence in protein_sequence:
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
            writer.writerow([sequence, bitstring, len(bitstring), mRNA_sequence, round(value.real, 3), executing_time])

            index += 1

def count_gates(protein_sequence):
    index = 1
    count = len(protein_sequence)

    results = []
    results.append(["amino_sequence", "qubit_len", "Gate_counts"])
    for sequence in protein_sequence:
        tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
        codonOpt = CodonOptimizationOneHot(sequence, tune_parameters)
        qubit_op = codonOpt.create_qubit_op()
        print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
        results.append([sequence, codonOpt.qubit_len, util.count_pauli_z(qubit_op)])
        index += 1
    return results

protein_sequence = "KIILFLALITLATCELYHYQECVRGTTVLLKEPCSSGTYEGNSPFHPLADNKFALTCFSTQFAFACPDGVKHVYQLRARSVSPKLFIRQEEVQELYSPIFLIVAAIVFITLCFTLKRKTE"
protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
for i in range(19, 20):
    split_len = i
    encoding_type = "one_hot"
    sequence_list = util.spilt_sequence(protein_sequence, split_len)
    name = time.strftime("%m-%d_%H-%M-%S")
    # write_data(f"../output/{split_len}/{encoding_type}/P0DTC2_{name}.csv", sequence_list[10:])
    data = count_gates(sequence_list)
    util.write_data(f"../output/{split_len}/Gate/{encoding_type}_P0DTC2_Gate.csv", data)


