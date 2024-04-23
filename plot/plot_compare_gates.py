import matplotlib.pyplot as plt
import time, csv


def plot_gates(dense_gates, one_hot_gates, output_file):
    x = range(6, 20)
    plt.plot(x, one_hot_gates[3:], marker='o', label='one-hot')
    plt.plot(x, dense_gates[3:], marker='*', label='dense')
    plt.xlabel("the length of each fragment of protein")
    plt.ylabel("the required number of gates")
    plt.title("the number required of gates with one-hot and dense encoding")
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

dense_gates, one_gates = read_data(f"../output/statistic_gates.csv")
plot_gates(dense_gates, one_gates, f"../output/figure_gates.png")


