import matplotlib.pyplot as plt
import time, csv


def plot_qubits(dense_qubits, one_hot_qubits, output_file):
    x = range(6, 20)
    plt.plot(x, one_hot_qubits[0], marker='o', label='maximum on one-hot')
    plt.plot(x, dense_qubits[0], marker='*', label='maximum on dense')
    plt.plot(x, one_hot_qubits[1], marker='h', label='mean on one-hot')
    plt.plot(x, dense_qubits[1], marker="v", label="mean on dense")
    plt.xlabel("the length of each fragment of protein")
    plt.ylabel("the number of qubits")
    plt.title("the required number of gates with one-hot and dense encoding")
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
dense = [dense_qubits, mean_dense_qubits]
one_hot = [one_hot_qubits, mean_one_hot_qubits]
plot_qubits(dense, one_hot, f"../output/figure_qubits.png")


