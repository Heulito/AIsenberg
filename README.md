# AIsenberg

# Bioactivity prediction to inhibit ERK2 enzyme
Prediction of molecules (*SMILES*) towards the inhibitation of ERK2 enzyme

# Start the web Application
1. Clone Project
```
git clone git@github.com:Heulito/AIsenberg.git
```
2. cd into project
```
cd AIsenberg
```
3. Create conda environment called WebAPP
```
conda create -n WebAPP python=3.8
```
4. Activate conda environment
```
conda activate WebAPP
```
5. Install the dependencies
```
pip install -r requirements.txt
```
7. run the Web Application
```
streamlit run WebAPP_IC50.py
```
