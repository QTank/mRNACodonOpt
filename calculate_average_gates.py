import os, csv
import util
import numpy as np


def get_average_gates(file_name):
    data = []
    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar="|")
        for row in reader:
            if row['Gate_counts'] == None:
                continue
            data.append(int(row["Gate_counts"]))

    return int(np.ceil(np.mean(data)))



path = "output"

data = []
data.append(["qubit_len", "dense", "one_hot"])
for i in range(3, 20):
    dense_average_gates, one_hot_average_gates = 0, 0
    file_list = os.listdir(f"{path}/{i}/Gate")
    for file_name in file_list:
        if "dense" in file_name:
            dense_average_gates = get_average_gates(f"{path}/{i}/Gate/{file_name}")

        if "one_hot" in file_name:
            one_hot_average_gates = get_average_gates(f"{path}/{i}/Gate/{file_name}")

    data.append([i, dense_average_gates, one_hot_average_gates])

util.write_data(f"{path}/statistic_gates.csv", data)



