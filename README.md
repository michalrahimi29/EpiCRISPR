# EpiCRISPR
We share here EpiCRISPR which is a tool used to predict CRISPR/Cas9 on-target efficiency by adding epigenetic and high-throughput data information. EpiCRISPR is a multi-layered preceptron, i.e. a cascade of fully connected layers of 64, 8 neurons, terminating with a final neuron, which outputs the on-target efficiency. The input is a fixed RNA sequence of 32nt length, ten possible epigenetic marks: Chromatin accessibility, CTCF binding, DNA methylation and seven histone modification: H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K36me3, H3K9ac and H3K27ac, and information from high-throughput data based model (We used CRISPROn scores prediction).
We utilized the Leenay et al. dataset to evaluate the contribution of flanking sequences epigenetic marks and high-throughput data information to on-target efficiency prediction and found that epigenetic marks play crucial role in predicting the on-target efficiency of CRISPR/Cas9. Moreover, we successfully show the generalizability of our trained model by successfully predicting the on-target efficiency to other cell types. 

# Requirements
The model is implemented with Keras 2.8.4 using Tensorflow backend and numpy 1.24.3.

# Run
To run EpiCRISPR follow the instruction bellow:

### Model Training and evaluation:
For traning the model and evaluating it with specific epigenetic marks utilizing the Leenay et al. dataset: First input should be 1 and afterward give the epigenetic marks you would like to use: chromatin_accessibility, CTCF_binding, H3K4me3, methylation. 


For training and evaluating the model on the Leenay et al. dataset without epigentic marks:
```python
python run.py 1
```
For training and evaluating the model on the Leenay et al. dataset with the epigentic marks- chromatin_accessibility, CTCF_binding:
```python
python run.py 1 chromatin_accessibility CTCF_binding
```
### Model evaluation on other human cell types:
For evaluating the model on the human cell types: HEK293, HCT116, K562 the first input should be the wanted cell type and and afterward give the epigenetic marks you would like to use: chromatin_accessibility, CTCF_binding, H3K4me3, methylation.


For training the model on the Leenay et al. dataset and evaluating on K562 with all epigenetic marks:
```python
python run.py K562 chromatin_accessibility CTCF_binding H3K4me3 methylation
```
For training the model on the Leenay et al. dataset and evaluating on HCT116 with epigenetic mark: chromatin_accessibility:
```python
python run.py HCT116 chromatin_accessibility
```
### Model for further use:
For further use with the trained model on the entire Leenay et al. dataset with the four epigenetic marks use the following code:
```python
reconstructed_model = keras.models.load_model("CRISPRepi_model.keras")
```
