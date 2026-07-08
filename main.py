gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

def read_fasta(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()

    genome = ""

    for line in lines:
        line = line.rstrip("\n")

        # Skip FASTA header line
        if line.startswith(">"):
            continue

        genome = genome + line.strip()

    return genome.upper()


def translate_genome(genome, frame):
    protein = ""

    # frame will be 0, 1, or 2
    for i in range(frame, len(genome) - 2, 3):
        codon = genome[i:i+3]
        amino_acid = gencode.get(codon, "X")
        protein = protein + amino_acid

    return protein


def find_long_proteins(translated_protein, frame, min_length=100):
    results = []

    pieces = translated_protein.split("*")

    aa_position = 0

    for piece in pieces:
        length = len(piece)

        if length >= min_length:
            # Convert amino acid position back to DNA position
            dna_start = frame + (aa_position * 3) + 1
            dna_end = frame + ((aa_position + length) * 3)

            results.append({
                "frame": frame + 1,
                "start": dna_start,
                "end": dna_end,
                "length": length,
                "protein": piece
            })

        # +1 because split removed the stop codon
        aa_position = aa_position + length + 1

    return results


# Change this to your actual file path
filename = "/share/SARS-class/SARS-2020.fasta"

genome = read_fasta(filename)

all_proteins = []

for frame in range(3):
    translated = translate_genome(genome, frame)
    long_proteins = find_long_proteins(translated, frame, 100)
    all_proteins = all_proteins + long_proteins

print("Genome length:", len(genome))
print("Number of long proteins found:", len(all_proteins))
print()

for i in range(len(all_proteins)):
    protein = all_proteins[i]

    print("Protein", i + 1)
    print("Frame:", protein["frame"])
    print("Start:", protein["start"])
    print("End:", protein["end"])
    print("Length:", protein["length"])
    print("Sequence:", protein["protein"])
    print()


# Save proteins to a FASTA file for BLAST
out = open("proteins_found.fasta", "w")

for i in range(len(all_proteins)):
    protein = all_proteins[i]

    out.write(">protein_" + str(i + 1))
    out.write("_frame_" + str(protein["frame"]))
    out.write("_start_" + str(protein["start"]))
    out.write("_end_" + str(protein["end"]))
    out.write("_length_" + str(protein["length"]))
    out.write("\n")

    out.write(protein["protein"])
    out.write("\n")

out.close()

print("Saved proteins to proteins_found.fasta")
