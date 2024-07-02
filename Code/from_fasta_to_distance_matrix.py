from Bio import SeqIO
import numpy as np
import csv

def hamming_distance(seq1, seq2):
    """计算两个序列之间的汉明距离"""
    if len(seq1) != len(seq2):
        raise ValueError("sequence length not match")
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def calculate_distance_matrix(fasta_file):
    """从FASTA文件读取序列，并计算汉明距离矩阵"""
    sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            if i != j:
                distance_matrix[i, j] = hamming_distance(sequences[i], sequences[j])

    return distance_matrix

# please change this fasta_file to the real path
print('path of fasta file required:')
fasta_file = input()
distance_matrix = calculate_distance_matrix(fasta_file)
print(distance_matrix)

def save_matrix_to_csv(matrix, file_path):
    # Open the file in write mode
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the number of rows as the first row
        writer.writerow([matrix.shape[0]])

        # Write the matrix rows
        writer.writerows(matrix)

# please change this output_csv_path to the real path
print('path of output file required:')
output_csv_path = input()

# Save the filled matrix to a CSV file
save_matrix_to_csv(distance_matrix, output_csv_path)
