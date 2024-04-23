from exact_bit import ExactBit, get_opt_bitstring
import util, time
from tunable_parameters import TunableParameters

def calculate_one_amino(amino_sequence):
    start_time = time.time()
    tunable_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
    opt_bitstring, opt_value = get_opt_bitstring(amino_sequence, tunable_parameters)
    mRNA_sequence = util.decode_bitstring(amino_sequence, opt_bitstring)
    executing_time = round(time.time() - start_time, 2)
    return [amino_sequence, opt_bitstring, len(opt_bitstring), mRNA_sequence, opt_value, executing_time]

def get_all(protein_sequence, split_len):

    results = []
    sequence_list = util.spilt_sequence(protein_sequence, split_len)
    count = len(sequence_list)
    for index in range(count):
        amino_sequence = sequence_list[index]
        data = calculate_one_amino(amino_sequence)
        print(f"{index}/{count} Start to calculate {amino_sequence} ... , requiring {util.required_qubit_len(amino_sequence)} qubits")
        results.append(data)
    return results


def write_exact_all(protein_sequence, split_len):
    path = f"../output/{split_len}/exact"
    name = time.strftime("%m-%d_%H-%M-%S")

    name = f"{path}/P0DTC2_{name}.csv"
    results = get_all(protein_sequence, split_len)
    util.write_data(name, results)


if __name__ == '__main__':

    protein_sequence = "FVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
    split_len = 3
    write_exact_all(protein_sequence, split_len)




