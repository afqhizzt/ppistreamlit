import streamlit as st
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
from Bio import SeqIO # I tambah ni untuk protein length & weight line 45-51
import json
import random  

def retrieve_ppi(protein_identifier, is_uniprot=True):
    if is_uniprot:
        string_url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": protein_identifier,
            "species": 9606  
        }
    else: # tak dapat
        string_url = "https://string-db.org/api/json/sequence"
        params = {
            "sequence": protein_identifier,
            "species": 9606  
        }
    response = requests.get(string_url, params=params)
    data = response.json()  
    network_df = pd.json_normalize(data)
    return network_df


def visualize_ppi(network_df):
    fig, ax = plt.subplots()
    ppi_graph = nx.from_pandas_edgelist(network_df, "preferredName_A", "preferredName_B")
    
    node_colors = [random.choice(['red', 'blue', 'green', 'yellow', 'orange', 'pink', 'cyan', 'purple']) for _ in range(len(ppi_graph.nodes))]
    graph = nx.spring_layout(ppi_graph)
    
    nx.draw(ppi_graph, graph, with_labels=True, node_size=500, node_color=node_colors, font_size=8, ax=ax)
    st.pyplot(fig) # sama macam plt.show()


def display_characteristics(network_df):
    ppi_graph = nx.from_pandas_edgelist(network_df, "preferredName_A", "preferredName_B")

    # Calculating protein length and weight
    protein_sequence = network_df["preferredName_A"][0]  # Taking the sequence of the first protein as an example
    protein_length = len(protein_sequence)

    # Check if the sequence is a valid protein sequence
    if all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in protein_sequence):
        protein_weight = molecular_weight(Seq(protein_sequence), "protein")
    else:
        protein_weight = "Sequence is not a valid protein sequence"

    st.write("Protein Length:", protein_length)
    st.write("Protein Weight:", protein_weight)

    st.write("Number of edges:", ppi_graph.number_of_edges())  # number of interactions
    st.write("Number of nodes:", ppi_graph.number_of_nodes())  # number of proteins
    st.write("Number of interactions for each node:", ppi_graph.degree())  # number of interactions for each node
    
def display_centrality(network_df):
    ppi_graph = nx.from_pandas_edgelist(network_df, "preferredName_A", "preferredName_B")
    
    # degree centrality
    degree_central = nx.degree_centrality(ppi_graph)
    st.write("Degree Centrality: ", degree_central)

    # top 5 proteins with highest centrality
    top_5_proteins = sorted(degree_central.items(), key=lambda x:-x[1])[:5]
    st.write("Top 5 Proteins with Highest Degree Centrality:")
    for protein, centrality in top_5_proteins:
        st.write(f"{protein}: {centrality}")

    high_centrality_proteins = [protein for protein, _ in top_5_proteins]
    
    # colors for the nodes in the main graph
    main_node_colors = [random.choice(['red', 'blue', 'green', 'yellow', 'orange', 'pink', 'cyan', 'purple']) for _ in range(len(ppi_graph.nodes))]
    # highlight color that is different from the main node colors
    highlight_color = random.choice(['red', 'blue', 'green', 'yellow', 'orange', 'pink', 'cyan', 'purple'])
    while highlight_color in main_node_colors:
        highlight_color = random.choice(['red', 'blue', 'green', 'yellow', 'orange', 'pink', 'cyan', 'purple'])

    fixed_layout = nx.spring_layout(ppi_graph)

    fig, ax = plt.subplots()
    nx.draw(ppi_graph, fixed_layout, with_labels=True, node_size=500, node_color=main_node_colors, font_size=8, ax=ax)
    nx.draw_networkx_nodes(ppi_graph, fixed_layout, nodelist=high_centrality_proteins, node_color=highlight_color, ax=ax)
    
    st.pyplot(fig)


def main():
    st.title('Protein Data Analysis App')

    st.write('**Instructions:**')
    st.write('1. Select either Uniprot ID or a protein sequence from the dropdown.')
    st.write('2. Enter a Uniprot ID or a protein sequence.')
    st.write('3. Press the "Retrieve Data" button to fetch protein-protein interaction network, its characteristics and centrality.')

    input_option = st.selectbox('Select Input Option', ['Uniprot ID', 'Protein Sequence'])
    protein_identifier = st.text_input('Enter Uniprot ID or Protein Sequence')

    if st.button('Retrieve Data'):
        if input_option == 'Uniprot ID':
            network_df = retrieve_ppi(protein_identifier)
        else:
            network_df = retrieve_ppi(protein_identifier, is_uniprot=False)

        # Tabs for displaying Protein-Protein Interaction Network, Protein Characteristics, and Centrality
        tabs = st.tabs(["Protein-Protein Interaction Network", "Protein Characteristics", "Protein Centrality"])

        with tabs[0]:
            st.subheader('Protein-Protein Interaction Network')
            visualize_ppi(network_df)

        with tabs[1]:
            st.subheader('Protein Characteristics')
            display_characteristics(network_df)

        with tabs[2]:
            st.subheader('Protein Centrality')
            display_centrality(network_df)

if __name__ == '__main__':
    main()
