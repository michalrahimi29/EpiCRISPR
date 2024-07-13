# EpiCRISPR
We share here EpiCRISPR which is a tool used to predict CRISPR/Cas9 on-target efficiency by adding epigenetic and high-throughput data information. EpiCRISPR is a multi-layered preceptron, i.e. a cascade of fully connected layers of 64, 8 neurons, terminating with a final neuron, which outputs the on-target efficiency. The input is a fixed RNA sequence of 32nt length, ten possible epigenetic marks: Chromatin accessibility, CTCF binding, DNA methylation and seven histone modification: H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K36me3, H3K9ac and H3K27ac, and information from a high-throughput data-based model (we used CRISPROn's score predictions).
We utilized the Leenay et al. dataset to evaluate the contribution of flanking sequences, epigenetic marks and high-throughput data information to on-target efficiency prediction and found that epigenetic marks play crucial role in predicting the on-target efficiency of CRISPR/Cas9. Moreover, we successfully show the generalizability of our trained model by successfully predicting the on-target efficiency to other cell types. 

# Requirements
The model is implemented with Keras 2.8.4 using Tensorflow backend and numpy 1.24.3.

# Run
To run EpiCRISPR follow the instruction bellow:

### Model Training and evaluation:
For traning the model and evaluating it with specific epigenetic markers and high-throughput data information, utilizing the Leenay et al. dataset: First input should be 1, afterward for adding CRISPROn's score predictions the second input shold be True else False, and last give the epigenetic markers you would like to use: chromatin, CTCF, methylation, H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K36me3, H3K9ac and H3K27ac. If you desire to use all epigenetic markers the third input should be 'all'.


For training and evaluating the model on the Leenay et al. dataset without CRISPROn scores and with chromatin accessibility marker:
```python
python run.py 1 False chromatin
```
For training and evaluating the model on the Leenay et al. dataset with CRISPROn scores and with CTCF binding and DNA methylation markers:
```python
python run.py 1 True CTCF methylation
```
For training and evaluating the model on the Leenay et al. dataset with CRISPROn score and all epigenetic markers:
```python
python run.py 1 True all
```
### Model evaluation on other human cell types:
For evaluating the model on the human cell types- HEK293, HCT116, K562 and T cell: First input should be 2, afterward for adding CRISPROn's score predictions the second input shold be True else False, the third input should be the wanted cell type and last give the epigenetic markers you would like to use, if you desire to use all avalible epigenetic markers the third input should be 'all'. For evaluation on other human cell type only the following epigenetic markers are available: chromatin, CTCF, H3K4me1, H3K27ac, H3K36me3. 


For training the model on the Leenay et al. dataset and evaluating on K562 without CRISPROn scores and with all epigenetic markers:
```python
python run.py 2 K562 False all
```
For training the model on the Leenay et al. dataset and evaluating on T-cell with CRISPROn scores and with H3K9me3 chromatin accessibility and CTCF binding markers:
```python
python run.py T True H3K9me3 chromatin CTCF
```
### Model for further use:
For further use with the trained model on the entire Leenay et al. dataset with all epigenetic markers and CRISPROn scores use the following code:
```python
reconstructed_model = keras.models.load_model("EpiCRISPR_model.keras")
```
