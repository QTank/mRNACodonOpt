import python_codon_tables as pct
import numpy as np
import util


class CodonOneHot:

    def __init__(self, amino, start_index, qubit_len, encoding_type=0, table_name="e_coli_316407"):
        """
        :param amino: the name of amino
        :param start_index: the position in the protein sequence
        :param qubit_len: the required number of qubits
        :param encoding_type: the type of encoding codon, 0 represents one-hot encoding, 1 represents dense encoding
        :param table_name: the name of codon table
        """

        self.amino = amino
        self.qubit_len = qubit_len
        self.encoding_type = encoding_type
        self.start_index = start_index
        self.codon_table = pct.get_codons_table(table_name)
        self.encoding_qubit_len = len(self.codon_table[amino])
        self.indicator_list = self._create_codon_pauli_z()

        self.indicator_dict = self.indicator()
        self.redundant_indicator_list = self.create_redundant_indicator()

    def _create_codon_pauli_z(self):
        return [util.build_indicator_qubit(self.qubit_len, self.start_index+i) for i in range(self.encoding_qubit_len)]

    def indicator(self):
        indicator_dict = {}

        for index, codon in enumerate(self.codon_table[self.amino]):
            indicator_dict[codon] = (self.indicator_list[index], self.codon_table[self.amino][codon])
        return indicator_dict

    def create_redundant_indicator(self):
        redundant_indicator_list = []
        codon_index = [2 ** i for i in range(self.encoding_qubit_len)]

        for index in range(0, 2 ** self.encoding_qubit_len):
            if index not in codon_index:
                redundant_indicator_list.append(self._build_indicator(index))

        return redundant_indicator_list

    def _build_indicator(self, index):
        binary_str = util.int_to_binary(index, self.encoding_qubit_len)
        indicator_qubit = util.build_full_identity(self.qubit_len)

        for i in range(len(binary_str)):
            if binary_str[i] == "0":
                indicator_qubit = indicator_qubit @ (util.build_full_identity(self.qubit_len) - self.indicator_list[i])
            else:
                indicator_qubit = indicator_qubit @ self.indicator_list[i]

        return indicator_qubit

