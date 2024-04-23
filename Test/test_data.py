import util

spilt_len = 3
dense = "dense"
exact = "exact"
one_hot = "one_hot"
dense_name = f"../output/{spilt_len}/P0DTC2_04-19_12-04-21.csv"
one_hot_name = f"../output/{spilt_len}/{one_hot}/P0DTC2_04-19_10-06-14.csv"
exact_name = f"../output/{spilt_len}/{exact}/P0DTC2_04-19_13-13-13.csv"
dense_dict = util.read_data(dense_name)
exact_dict = util.read_data(exact_name)
one_hot_dict = util.read_data(one_hot_name)


for sequence in dense_dict:

    if float(dense_dict[sequence]['value']) - float(one_hot_dict[sequence]['value']) > 0.1:
        print(sequence, dense_dict[sequence]['value'], one_hot_dict[sequence]['value'])
        print("dense and one hot mRNA_sequence equal: ", one_hot_dict[sequence]["mRNA_sequence"] == one_hot_dict[sequence]["mRNA_sequence"])
        print(f"exact:{exact_dict[sequence]}, dense:{dense_dict[sequence]}, one hot:{one_hot_dict[sequence]}")