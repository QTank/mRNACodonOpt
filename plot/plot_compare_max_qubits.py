import matplotlib.pyplot as plt
import time, csv


def plot_qubits(dense_qubits, one_hot_qubits, output_file):
    x = range(6, 17)
    plt.plot(x, one_hot_qubits[0], marker='o', label='One-Hot Encoding')
    plt.plot(x, dense_qubits[0], marker='*', label='Dense Encoding')
    plt.xlabel("The length of each amino acid sequence")
    plt.ylabel("The required max number of qubits")
    plt.title("The required max number of qubits with one-hot and dense encoding")
    plt.legend()
    plt.savefig(output_file)


def read_data(file_name):
    dense_gates = []
    one_gates = []
    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar="|")
        for row in reader:
            dense_gates.append(int(row['dense']))
            one_gates.append(int(row['one_hot']))

    return dense_gates, one_gates

dense_qubits, one_hot_qubits = read_data(f"../output/statistic_qubits.csv")
mean_dense_qubits, mean_one_hot_qubits = read_data(f"../output/statistic_mean_qubits.csv")
dense = [dense_qubits[:-3], mean_dense_qubits[:-3]]
one_hot = [one_hot_qubits[:-3], mean_one_hot_qubits[:-3]]
plot_qubits(dense, one_hot, f"../output/max_figure_qubits.png")


