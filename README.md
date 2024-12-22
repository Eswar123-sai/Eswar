from Bio import Entrez, SeqIO
from Bio import BiopythonParserWarning  # Import BiopythonParserWarning explicitly
import warnings

# Suppress warnings from Biopython
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

# Set email for Entrez
Entrez.email = "deswar902@gmail.com"  # Required by NCBI

# Function to fetch sequences based on taxonomy ID and database type

def fetch_sequences(tax_id, db_type, retmax=100):
    try:
        handle = Entrez.esearch(db=db_type, term=f"txid{tax_id}[Organism]", retmax=retmax)
        record = Entrez.read(handle)
        ids = record["IdList"]
        
        sequences = []
        for seq_id in ids:
            handle = Entrez.efetch(db=db_type, id=seq_id, rettype="gb", retmode="text")
            seq_record = SeqIO.read(handle, "genbank")
            sequences.append(seq_record)
        return sequences
    except Exception as e:
        print(f"Error fetching sequences: {e}")
        return []
 
# Fetch protein and nucleotide sequences for Solanum lycopersicum
Solanum_lycopersicum_tax_id = "4081"
protein_sequences = fetch_sequences(Solanum_lycopersicum_tax_id, "protein")
gene_sequences = fetch_sequences(Solanum_lycopersicum_tax_id, "nucleotide")

# Display the number of sequences fetched
print(f"Fetched {len(protein_sequences)} protein sequences.")
print(f"Fetched {len(gene_sequences)} gene sequences.")

# Save the sequences to FASTA files
SeqIO.write(protein_sequences, "Solanum_lycopersicum_proteins.fasta", "fasta")
SeqIO.write(gene_sequences, "Solanum_lycopersicum_genes.fasta", "fasta")
print("Sequences saved to FASTA files: Solanum_lycopersicum_proteins.fasta and Solanum_lycopersicum_genes.fasta")

def remove_redundancy(sequences):
    unique_sequences = []
    seen_sequences = set()
    
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        if seq_str not in seen_sequences:
            seen_sequences.add(seq_str)
            unique_sequences.append(seq_record)
    
    return unique_sequences

# Remove redundancy from protein and gene sequences
unique_protein_sequences = remove_redundancy(protein_sequences)

# Save the unique sequences to new FASTA files
SeqIO.write(unique_protein_sequences, "unique_Solanum_lycopersicum_proteins.fasta", "fasta")

print("Unique sequences saved to FASTA file: unique_Solanum_lycopersicum_proteins.fasta")

for rec in gene_sequences:
    fets = rec.features
    print(len(fets))

    import matplotlib.pyplot as plt

def construct_gene_structure(gene_record):
    features = []
    for feature in gene_record.features:
        if feature.type == "CDS":
            features.append({"start": feature.location.start, "end": feature.location.end, "type": "exon"})
        else:
            features.append({"start": feature.location.start, "end": feature.location.end, "type": feature.type})
    return features

# Visualize the gene structure
def visualize_gene_structure(gene_structure):
    plt.figure(figsize=(10, 2))
    for feature in gene_structure:
        plt.plot([feature["start"], feature["end"]], [1, 1], label=feature["type"], lw=6 if feature["type"] == "exon" else 2)
    plt.title("Gene Structure")
    plt.show()

# Example of constructing and visualizing gene structure for one gene
gene_structure = construct_gene_structure(gene_sequences[0])
visualize_gene_structure(gene_structure)

import matplotlib.pyplot as plt

# Function to construct gene structure
def construct_gene_structure(gene_record):
    features = []
    for feature in gene_record.features:
        if feature.type == "CDS":  # Focus on the CDS (coding sequences)
            features.append({
                "start": int(feature.location.start),  # Extract integer start position
                "end": int(feature.location.end),      # Extract integer end position
                "type": "exon"
            })
        else:
            features.append({
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "type": feature.type
            })
    return features

# Function to visualize gene structure
def visualize_gene_structure(gene_structure, gene_id):
    plt.figure(figsize=(10, 2))
    for feature in gene_structure:
        plt.plot([feature["start"], feature["end"]], [1, 1], label=feature["type"], lw=6 if feature["type"] == "exon" else 2)
    plt.title(f"Gene Structure for Gene ID: {gene_id}")
    plt.show()

# Loop through all genes and visualize their structure
for gene in gene_sequences:
    gene_structure = construct_gene_structure(gene)  # Construct the gene structure
    visualize_gene_structure(gene_structure, gene.id)  # Visualize the structure for each gene

 import re

def search_motifs(sequence, motif_pattern):
    matches = re.finditer(motif_pattern, str(sequence))
    motifs = []
    for match in matches:
        motifs.append((match.start(), match.end(), match.group()))
    return motifs

# Define a sample motif pattern (e.g., serine/threonine-rich regions)
motif_pattern = r"[ST]P"

# Example usage of the search_motifs function with a list of sequences
for i, protein_seq in enumerate(unique_protein_sequences):
    motifs = search_motifs(protein_seq.seq, motif_pattern)
    print(f"Motifs found in Protein Sequence {i+1}: {motifs}")

with open("motif_search_results.txt", "w") as file:
    for i, protein_seq in enumerate(unique_protein_sequences):
        motifs = search_motifs(protein_seq.seq, motif_pattern)
        file.write(f"Motifs found in Protein Sequence {i+1}:\n")
        for start, end, match in motifs:
            file.write(f" - Start: {start}, End: {end}, Sequence: {match}\n")
        file.write("\n")  # Add a newline for readability

print("Motif search results saved to 'motif_search_results.txt'")

import os
import subprocess
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from Bio import SeqIO

# Load protein sequences
fasta_file = "Solanum_lycopersicum_proteins.fasta"  # Make sure this file exists and is correctly specified
setaria_proteins = list(SeqIO.parse(fasta_file, "fasta"))

# File for alignment output
alignment_file = "aligned_sequences.aln"

# Run MUSCLE alignment
try:
    result = subprocess.run(
        ["C:/Program Files/Tools/Muscle/muscle-windows-v5.2.exe", "-align", fasta_file, "-output", alignment_file],
        check=True,
        text=True,
        capture_output=True
    )
    print("MUSCLE output:", result.stdout)
except subprocess.CalledProcessError as e:
    print("Error running MUSCLE:", e.stderr)
    raise

# Step 2: Read the alignment
alignment = AlignIO.read(alignment_file, "fasta")

# Step 3: Calculate the distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Step 4: Construct the phylogenetic tree using Neighbor-Joining
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

print(tree)

# Step 5: Visualize the tree
plt.figure(figsize=(10, 5))  # Set a larger figure size
Phylo.draw(tree, do_show=True)  # Draw the tree without showing it yet

# Optionally, adjust the limits of the axes to ensure all clades are visible
plt.xlim(0, 10)  # Adjust the x-axis limits as needed
plt.ylim(0, 20)  # Adjust the y-axis limits as needed

# Save the figure as a PDF
plt.savefig("phylogenetic_tree.pdf", format='pdf', bbox_inches='tight')  # Save as PDF
plt.close()  # Close the figure to free up memory

# Save the tree to a file
output_image_path = "phylogenetic_tree.jpg"
plt.savefig(output_image_path, dpi=300, bbox_inches="tight")  # Save as high-resolution PNG
plt.show()


print(f"Tree visualization saved to: {output_image_path}")

plt.figure(figsize=(10, 5))  # Set a larger figure size
Phylo.draw(tree, do_show=True)  # Draw the tree without showing it yet

plt.savefig(output_image_path, dpi=300, bbox_inches="tight")  # Save as high-resolution PNG
    
