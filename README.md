# mRNA codon optimization
This method is from the paper https://arxiv.org/pdf/2404.14858. For a detailed explanation, you can refer to the paper. 


## 1. Create virtual environment
The virtual environment can isolate dependencies and prevent version conflicts of dependencies between different projects. It ensures that the installed packages and libraries are exclusively used by this project, avoiding effects on other projects.

Refer to the following links about how to create virtual environment:
1. https://docs.python.org/3/library/venv.html
2. https://www.arch.jhu.edu/python-virtual-environments/

## 2. Download the code of mRNA codon optimization
Clone the mRNA codon optimization from Github:

git clone https://github.com/QTank/mRNACodonOpt.git

Note: make sure Git is installed git before cloning the code.

## 3. Install required packages
pip install -r requirement.txt

Note: make sure that pip is installed before running this command.

## 4. Run mRNA codon optimization
Navigate to the solution folder where the optimization files are located:

1. Dense Encoding way (solution.py)

   Use dense encoding way to optimize mRNA codons on quantum computing simulation in Qiskit
2. One-Hot Encoding way (solution_one_hot.py)
   Uses one-hot encoding way to optimize mRNA codons on quantum computing simulation in Qiskit
3. Exact Value optimization ('solution_exact.py')
   Uses brute force to compute the exact value of the Hamiltonian for mRNA codon optimization 


