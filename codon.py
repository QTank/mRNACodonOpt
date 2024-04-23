import python_codon_tables as pct
import numpy as np
import util


class Codon:

    def __init__(self, amino, start_index, qubit_len, encoding_type=1, table_name="e_coli_316407"):
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
        self.encoding_qubit_len = self._calculate_used_qubit()
        if self.encoding_qubit_len == 0:
            self.codon_qubit_list = []
            self.indicator_dict = self.only_one_codon()
            self.redundant_indicator_list = []
        else:
            self.codon_qubit_list = self._create_codon_pauli_z()
            self.indicator_dict, self.redundant_indicator_list = self.indicator()

    def _create_codon_pauli_z(self):
        return [util.build_indicator_qubit(self.qubit_len, self.start_index+i) for i in range(self.encoding_qubit_len)]

    def indicator(self):
        indicator_dict = {}
        redundant_indicator_list = []

        for index, codon in enumerate(self.codon_table[self.amino]):
            indicator_dict[codon] = (self._build_indicator(index), self.codon_table[self.amino][codon])

        for index in range(len(self.codon_table[self.amino]), 2 ** self.encoding_qubit_len):
            redundant_indicator_list.append(self._build_indicator(index))
        return indicator_dict, redundant_indicator_list

    def _calculate_used_qubit(self):
        if self.encoding_type == 0:
            return len(self.codon_table[self.amino])

        return round(np.log2(len(self.codon_table[self.amino])))

    def _build_indicator(self, index):
        binary_str = util.int_to_binary(index, self.encoding_qubit_len)
        indicator_qubit = util.build_full_identity(self.qubit_len)

        for i in range(len(binary_str)):
            if binary_str[i] == "0":
                indicator_qubit = indicator_qubit @ (util.build_full_identity(self.qubit_len) - self.codon_qubit_list[i])
            else:
                indicator_qubit = indicator_qubit @ self.codon_qubit_list[i]

        return indicator_qubit

    def only_one_codon(self):
        indicator_dict = {}
        for codon in self.codon_table[self.amino]:
            indicator_dict[codon] = (util.build_full_identity(self.qubit_len), self.codon_table[self.amino][codon])
        return indicator_dict
