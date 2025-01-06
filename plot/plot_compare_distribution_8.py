import matplotlib.pyplot as plt
from collections import Counter
import util


def plot_qubits(dense_dict, one_hot_dict, output_file):
    one_hot_qubits = list(one_hot_dict.keys())
    dense_qubits = list(dense_dict.keys())
    min_qubits = min(dense_qubits)
    max_qubits = max(one_hot_qubits)
    x = range(min_qubits, max_qubits+1)
    one_hot_sequence = []
    dense_sequence = []
    for i in x:
        if i not in one_hot_dict:
            one_hot_sequence.append(0)
        else:
            one_hot_sequence.append(one_hot_dict[i])
        if i not in dense_dict:
            dense_sequence.append(0)
        else:
            dense_sequence.append(dense_dict[i])


    plt.bar(x, one_hot_sequence, label='one-hot')
    plt.bar(x, dense_sequence, bottom=one_hot_sequence, label='dense')
    plt.xlabel("The number of qubits")
    plt.ylabel("The number of sequences")
    plt.title("The number of sequences for each fragment length of 8")
    plt.legend()
    plt.savefig(output_file)


protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
split_len = 8
sequence_list = util.spilt_sequence(protein_sequence, split_len)
dense_qubits = util.count_qubits(sequence_list, "dense")
one_hot_qubits = util.count_qubits(sequence_list, "one_hot")
dense_dict = Counter(dense_qubits)
one_hot_dict = Counter(one_hot_qubits)
plot_qubits(dense_dict, one_hot_dict, f"../output/qubits_distribution.png")

