import io
import streamlit as st
import requests
from Bio import SeqIO, SeqUtils
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Function to fetch protein-protein interaction network from STRING DB
def fetch_ppi_network(uniprot_id):
    url = f"https://string-db.org/api/json/interaction_partners?identifiers={uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for 4xx and 5xx status codes
        st.write("Response content:", response.content)
        data = response.json()
        st.write("Data structure:", data)
        interaction_partners = [partner["preferredName"] for partner in data if "preferredName" in partner]
        return interaction_partners
    except Exception as e:
        st.error("An unexpected error occurred while fetching protein-protein interaction network:", e)
        st.error("Response status code:", response.status_code)
        st.error("Response content:", response.content)
        st.error("Response text:", response.text)
        return None


# Function to fetch protein-protein interaction network from STRING DB
def fetch_ppi_network(uniprot_id):
    url = f"https://string-db.org/api/json/interaction_partners?identifiers={uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for 4xx and 5xx status codes
        data = response.json()
        interaction_partners = [partner["preferredName"] for partner in data[uniprot_id]]
        return interaction_partners
    except Exception as e:
        st.error("An unexpected error occurred while fetching protein-protein interaction network:", e)
        st.error("Response status code:", response.status_code)
        st.error("Response content:", response.content)
        st.error("Response text:", response.text)
        return None


# Function to perform sequence alignment
def perform_sequence_alignment(protein_sequence):
    # Write the protein sequence to a temporary FASTA file
    with open("temp.fasta", "w") as f:
        f.write(">query\n")
        f.write(protein_sequence)

    # Perform sequence alignment using Clustal Omega
    cline = ClustalOmegaCommandline(infile="temp.fasta", outfile="alignment.fasta", verbose=True, auto=True)
    stdout, stderr = cline()

    # Read the alignment output
    alignment = []
    with open("alignment.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            alignment.append(record.seq)
    
    return alignment

# Main function to run the Streamlit app
def main():
    st.title("Protein Data Analysis")

    option = st.sidebar.selectbox(
        'Choose Input Type',
        ('UniProt ID', 'Protein Sequence'))

    if option == 'UniProt ID':
        uniprot_id = st.sidebar.text_input('Enter UniProt ID')

        if st.sidebar.button('Fetch Data'):
            protein_data = fetch_protein_data(uniprot_id)
            if protein_data:
                st.write("### Protein Characteristics")
                st.write(f"Length: {protein_data['length']}")
                st.write(f"Molecular Weight: {protein_data['molecular_weight']}")

                st.write("### Protein-Protein Interaction Network")
                ppi_network = fetch_ppi_network(uniprot_id)
                if ppi_network:
                    st.write(ppi_network)
                else:
                    st.write("Failed to fetch PPI network.")

    elif option == 'Protein Sequence':
        protein_sequence = st.sidebar.text_area('Enter Protein Sequence')

        if st.sidebar.button('Analyze Sequence'):
            st.write("### Protein Characteristics")
            length = len(protein_sequence)
            molecular_weight = SeqUtils.molecular_weight(protein_sequence)
            st.write(f"Length: {length}")
            st.write(f"Molecular Weight: {molecular_weight}")

            st.write("### Sequence Alignment")
            alignment_output = perform_sequence_alignment(protein_sequence)
            st.write(alignment_output)

            st.write("### Additional Merits")
            protein_analysis = ProteinAnalysis(protein_sequence)
            isoelectric_point = protein_analysis.isoelectric_point()
            st.write(f"Isoelectric Point: {isoelectric_point}")

            hydrophobicity = protein_analysis.protein_scale(window=9, edge="center")[0]
            st.write(f"Hydrophobicity: {hydrophobicity}")

main()
