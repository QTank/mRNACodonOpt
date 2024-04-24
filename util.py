import numpy as np
from typing import Set
from qiskit.opflow import PauliOp, I, Z
from itertools import groupby
import python_codon_tables as pct
import time, csv

def build_full_identity(num_qubits: int) -> PauliOp:
    full_identity = I
    for _ in range(1, num_qubits):
        full_identity = I ^ full_identity
    return full_identity


def build_pauli_z_op(num_qubits: int, pauli_z_indices: Set[int]) -> PauliOp:
    if 0 in pauli_z_indices:
        operator = Z
    else:
        operator = I
    for i in range(1, num_qubits):
        if i in pauli_z_indices:
            operator = operator ^ Z
        else:
            operator = operator ^ I

    return operator


def build_indicator_qubit(qubit_len: int, pauli_z_index: int) -> PauliOp:
    return 0.5 * build_full_identity(qubit_len) - 0.5 * build_pauli_z_op(qubit_len, {pauli_z_index})


def int_to_binary(index, length):
    binary_str = bin(index)[2:]
    bin_len = len(binary_str)

    if length > bin_len:
        binary_str = "%s%s" % ("".join(["0" for i in range(length - bin_len)]), binary_str)

    return binary_str


def get_gc_count(codon):
    count = 0
    for ch in codon:
        if ch in ["G", "C"]:
            count += 1
    return count


def get_neighbour_repetition(c1, c2):
    count_dups = [sum(1 for _ in group) for _, group in groupby(c1+c2)]
    return max(count_dups)**2 - 1


def decode_bitstring(protein_sequence, bitstring, table_name='e_coli_316407'):
    mRNA_sequence = ""
    start_index = 0
    codon_table = pct.get_codons_table(table_name)
    try:
        for amino in protein_sequence:
            codon_list = list(codon_table[amino].keys())
            encoding_len = round(np.log2(len(codon_list)))

            if encoding_len == 0:
                mRNA_sequence += codon_list[0]
            else:
                loc = int(bitstring[start_index:start_index + encoding_len], 2)
                start_index += encoding_len
                mRNA_sequence += codon_list[loc]
    except:
        print(f"the bitstring {bitstring} is invalid!")

    return mRNA_sequence.replace("T", "U")


def spilt_sequence(protein_sequence, spilt_len):
    sequence_list = []

    index = 0
    while index < len(protein_sequence):
        sequence_list.append(protein_sequence[index:index+spilt_len])
        index += spilt_len

    return sequence_list


def decode_one_hot_bitstring(protein_sequence, bitstring, table_name='e_coli_316407'):
    mRNA_sequence = ""
    start_index = 0
    codon_table = pct.get_codons_table(table_name)
    try:
        for amino in protein_sequence:
            codon_list = list(codon_table[amino].keys())
            encoding_len = len(codon_list)
            valid_index = [2 ** i for i in range(encoding_len)]

            loc = int(bitstring[start_index:start_index + encoding_len], 2)
            if loc not in valid_index:
                print(f"the bitstring {bitstring} is invalid!")
                return mRNA_sequence

            codon_index = bitstring[start_index:start_index + encoding_len].find("1")
            mRNA_sequence += codon_list[codon_index]
            start_index += encoding_len
            # mRNA_sequence += codon_list[bitstring.find("1")]
    except:
        print(f"the bitstring {bitstring} is invalid!")

    return mRNA_sequence.replace("T", "U")


def write_data(name, result):
    with open(name, 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(result)

def decode_bitstring_list(protein_sequence, bitstring, table_name='e_coli_316407'):
    mRNA_sequence = []
    start_index = 0
    codon_table = pct.get_codons_table(table_name)
    try:
        for amino in protein_sequence:
            codon_list = list(codon_table[amino].keys())
            encoding_len = round(np.log2(len(codon_list)))

            if encoding_len == 0:
                mRNA_sequence.append(codon_list[0])
            else:
                loc = int(bitstring[start_index:start_index + encoding_len], 2)
                start_index += encoding_len
                mRNA_sequence.append(codon_list[loc])
    except:
        print(f"the bitstring {bitstring} is invalid!")

    return mRNA_sequence


def read_data(name):
    data = {}
    with open(name, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar="|")
        for row in reader:
            data[row['amino_sequence']] = {"value":row['value'], "mRNA_sequence":row["mRNA_sequence"]}
    return data


def encode_codon_sequence(protein_sequence, mrna_sequence, encoding_type, table_name='e_coli_316407'):

    codon_table = pct.get_codons_table(table_name)
    encoding_string = []
    if len(protein_sequence) * 3 != len(mrna_sequence):
        raise Exception("The length of protein sequence does not match the length of mRNA sequence!")

    for index in range(len(protein_sequence)):
        amino = protein_sequence[index]
        loc = 0
        for codon in codon_table[amino]:
            if codon == mrna_sequence:
                encoding_string.append()


def required_qubit_len(amino_sequence, table_name='e_coli_316407'):
    codon_table = pct.get_codons_table(table_name)
    qubit_len = 0
    for amino in amino_sequence:
        qubit_len += round(np.log2(len(codon_table[amino])))
    return qubit_len


def count_pauli_z(operator):
    count = 0
    if isinstance(operator, PauliOp):
        table_z = operator.primitive.z
        for i in range(len(table_z)):
            if table_z[i]:
                count += 1
        return count

    for hamiltonian in operator:
        table_z = hamiltonian.primitive.paulis.z[0]
        for i in range(len(table_z)):
            if table_z[i]:
                count += 1

    return count


def count_qubits(protein_sequence, encoding_type):
    results = []
    index = 1
    count = len(protein_sequence)
    for sequence in protein_sequence:
        qubit_len = calculate_qubit_len(sequence, encoding_type)
        print(f"{index}/{count} Start to execute VQE on {sequence} with dense encoding ... , requiring {qubit_len} qubits")
        index += 1
        results.append(qubit_len)
    print(f"The maximum required number of qubits: {max(results)}")

    return results


def calculate_qubit_len(sequence, encoding_type, table_name='e_coli_316407'):
    codon_table = pct.get_codons_table(table_name)
    count = 0
    for amino in sequence:
        if encoding_type == "one_hot":
            count += len(codon_table[amino])
        else:
            # count += np.ceil(np.log2(len(codon_table[amino])))
            count += round(np.log2(len(codon_table[amino])))

            if round(np.log2(len(codon_table[amino]))) != np.ceil(np.log2(len(codon_table[amino]))):

                print(sequence)

    return count
