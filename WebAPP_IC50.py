######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
import py3Dmol
import subprocess
import os
import base64
import pickle
from stmol import showmol
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from sympy import true



######################
# Molecular descriptor calculator
######################
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

################
# to download results
#################
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href


    

######################
# Page Title
######################
st.markdown("<h1 style='text-align: center; color: blue;'>WEB APPLICATION for IC50 Prediction</h1>", unsafe_allow_html=True)
st.write("""
This Web Application will allow you to predict the half maximal inhibitory concentrarion (**IC50**) for the ERK2 Enzyme.
The simplified pathway (left) shows how the ERK Enzyme is part of the regulation of cell grows [*1*]. If the ERK enzyme will be inhibited, cancer growth can be reduced [*2*]. 
The 3D - Visualisation of the ERK2 enzyme (right) will give an impression of the complexity of this Enzyme [*3*]. For further investigation, you are able to move and zoom into the enzyme[*4*]. 
""")
image = Image.open('Pathway.jpg')
new_image = image.resize((390,400))

#2ERK Protein
#If we have a Inhibitor which we want to show, we can look it up in 
#the www.rcsb.org datebase and visualize it 
xyzview = py3Dmol.view(query='pdb:2ERK') 
xyzview.setStyle({'cartoon':{'color':'spectrum'}})


#use st.columns to devide the page in two columns 
col1, col2 = st.columns([4,4])
with col1:
    st.image(new_image, caption="simplified Pathway")
with col2:
    showmol(xyzview, height = 500,width=800)    




st.write("""
# IC50 Prediction of Small Molecule

""")


######################
# Input molecules 
######################



def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,"mol")
    xyzview.setStyle({"stick":{}})
    xyzview.setBackgroundColor("white")
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=800)




####################
# Model building
####################
# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('erk_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    Smiles = pd.Series(load_data[0], name="Smiles")
    df = pd.concat([molecule_name, prediction_output, Smiles], axis=1).sort_values(by=['pIC50'], ascending=False)
    df = df[df["molecule_name"].str.contains("<NA>")==False]
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)
    


# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/Heulito/AIsenberg/main/predict_SMILES.txt)
""")
load_data = pd.read_table(uploaded_file, sep=' ', header=None)
load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

if st.sidebar.button('Predict'):
    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)

else:
    st.info('Upload input data in the sidebar to start!')

st.write("""
# Molecule structure 
You can investigate the molecule further by moving or zooming into it [*5*].
""")
Data = load_data.iloc[:, 0]
  
option = st.selectbox(
    'Which Molecule do you want to investigate?', Data, args=True)

st.write('You selected:', option)

# Sidebar, visualisation of the best molecule
#st.header("Please fill in the **SMILES** of the molecule you want to have a look at." )
#compound_smiles=st.text_input("SMILES input here:","CCCC")
blk=makeblock(option)
render_mol(blk) 
   







######################
# References
######################
st.write("""
# References

This *Web* *Application* and the used *Machin Learning Model* was inspired from Chanin Nantasenamat also known as the Data Professor.
https://www.youtube.com/c/DataProfessor/playlists/Bioinformatics%20Project
https://github.com/dataprofessor


[1] De Simone, A., Evanitsky, M. N., Hayden, L., Cox, B. D., Wang, J., Tornini, V. A., ... & Di Talia, S. (2021). Control of osteoblast regeneration by a train of Erk activity waves. Nature, 590(7844), 129-133.

[2] Yap JL, Worlikar S, MacKerell AD Jr, Shapiro P, Fletcher S. Small-molecule inhibitors of the ERK signaling pathway: Towards novel anticancer therapeutics. ChemMedChem. 2011 Jan 3;6(1):38-48. doi: 10.1002/cmdc.201000354. PMID: 21110380; PMCID: PMC3477473.

[3] Canagarajah, B.J., Khokhlatchev, A., Cobb, M.H., Goldsmith, E.J. Activation mechanism of the MAP kinase ERK2 by dual phosphorylation. (1997) Cell 90: 859-869, DOI: 10.1016/s0092-8674(00)80351-7

[4] Nápoles-Duarte JM, Biswas A, Parker MI, Palomares-Baez JP, Chávez-Rojo MA and Rodríguez-Valdez LM (2022) Stmol: A component for building interactive molecular visualizations within streamlit web-applications. Front. Mol. Biosci. 9:990846. doi: 10.3389/fmolb.2022.990846

[5] Molecular visualization in Streamlit using RDKit and Py3DMol by José Manuel Nàpoles Durate (2021)
https://github.com/napoles-uach/Medium_Mol


""")