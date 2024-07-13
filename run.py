import os
from EpiCRISPR import *
import tensorflow as tf
import random
import sys

"This file runs and initialize CRISPRepi"


def initialize(seed):
    randomNumSeed = 123
    np.random.seed(randomNumSeed + seed)  # Numpy random number generator
    tf.random.set_seed(0 + seed)  # TensorFlow random number generator
    random.seed(0 + seed)  # Python's built-in random number generator
    a = pd.read_csv("Final_leenay_dataset.csv")
    seqs_protospacer = a["protospacer"].to_numpy()
    seqs_pam = a["PAM"].tolist()
    seqs_up = a["upstream"].tolist()
    seqs_down = a["downstream"].tolist()
    no_var = a["no_variant"].to_numpy()
    reads = a["total_reads"].to_numpy()
    Total = np.sum(reads)
    W = np.divide(reads, Total)
    epsilon = 0.2
    weights = W + epsilon
    l = np.divide(no_var, reads)
    labels = np.subtract(1.0, l)
    epigeneticDic = {}
    flag = sys.argv[1]
    if flag == '1':  # 5CV on leenay dataset
        crispron = sys.argv[2]  # if True adding high-throughput data information
        epigenetics = sys.argv[3:]
    else:  # prediction on human cells
        cell_type = sys.argv[2]
        crispron = sys.argv[3]
        epigenetics = sys.argv[4:]
    if epigenetics == ['all']:
        if flag == '1' or (flag == '2' and cell_type == 'T'):
            files = [f for f in os.listdir('epigenetics/T/') if os.path.isfile(os.path.join('epigenetics/T/', f))]
            for file in files:
                key = file.split('_')[1].split('.')[0]
                b = pd.read_csv(f'epigenetics/T/{file}')
                epigeneticDic[key] = b['epigenetics'].to_numpy()
        else:
            files = [f for f in os.listdir(f'epigenetics/{cell_type}/') if os.path.isfile(os.path.join(f'epigenetics/{cell_type}/', f))]
            for file in files:
                key = file.split('_')[1]
                b = pd.read_csv(f'epigenetics/{cell_type}/{file}')
                epigeneticDic[key] = b['epigenetics'].to_numpy()
    else:
        for key in epigenetics:
            b = pd.read_csv(f'epigenetics/T/epigenetics_{key}.csv')
            epigeneticDic[key] = b['epigenetics'].to_numpy()
    if flag == '1':
        return lennaysRun(seqs_protospacer, seqs_pam, seqs_up, seqs_down, weights, labels, crispron, epigeneticDic)
    else:
        return lennayPredicionOnHumanCells(seqs_protospacer, seqs_pam, seqs_up, seqs_down, labels, cell_type, weights,
                                           crispron, epigeneticDic)
    # save_trained_model(seqs_protospacer, seqs_pam, seqs_up, seqs_down, weights, labels, crispron, identifiers)


# This function trains the model on all the Leenay dataset and save it. If one wants to test the model on different
# dataset he can use the trained model CRISPRepi_model.keras
def save_trained_model(seqs_protospacer, seqs_pam, seqs_up, seqs_down, weights, labels, crispron, epigeneticDic):
    data = createTrainSet(UP, DOWN, seqs_protospacer, seqs_down, seqs_up, seqs_pam, crispron, epigeneticDic)
    modeli = leenay_Model()
    modeli.fit(data, labels, epochs=25, batch_size=16, verbose=1, sample_weight=weights)
    modeli.save("EpiCRISPR_model.keras")


if __name__ == '__main__':
    initialize(135)
