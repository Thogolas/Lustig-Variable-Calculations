
import os
import subprocess
import pandas as pd

from Bio.Blast import NCBIXML
from entropy_calc import calculate_e6, calculate_e20


# Read sequences
sequences = pd.read_csv(
    "/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/entropy_rewrite/sequences.csv"
)


# Directories for BLAST files
blast_dir = (
    "/Users/patrickmincher/Desktop/LUSTIG SPSS/"
    "Generate_E6_E20_VKBAT/entropy_rewrite/blast_results"
)

os.makedirs(blast_dir, exist_ok=True)

db_path = "/Users/patrickmincher/blastdb/uniref50"

all_results = []

# Load checkpoint if it exists
checkpoint_file = os.path.join(blast_dir, "checkpoint.csv")
completed = set()
if os.path.exists(checkpoint_file):
    checkpoint_df = pd.read_csv(checkpoint_file)
    all_results = checkpoint_df.to_dict('records')
    completed = set(checkpoint_df[['Protein', 'Chain']].apply(lambda r: f"{r['Protein']}_{r['Chain']}", axis=1))
    print(f"Resuming from checkpoint - {len(completed)} proteins already done")


# Process each protein sequence
for idx, row in sequences.iterrows():

    protein = row["Protein"]
    chain = row["Chain"]
    seq = row["Sequence"]

    # Skip if already completed
    if f"{protein}_{chain}" in completed:
        continue

    print(
        f"BLASTing {protein} chain {chain} "
        f"({idx + 1}/{len(sequences)})"
    )


    # Write sequence to temporary FASTA file
    fasta_file = os.path.join(
        blast_dir,
        f"{protein}_{chain}.fasta"
    )

    with open(fasta_file, "w") as f:
        f.write(f">{protein}_{chain}\n{seq}\n")


    # Run local BLAST
    xml_file = os.path.join(
        blast_dir,
        f"{protein}_{chain}.xml"
    )

    subprocess.run([
        "blastp",
        "-query", fasta_file,
        "-db", db_path,
        "-out", xml_file,
        "-outfmt", "5",
        "-num_alignments", "1000",
        "-num_threads", "8"
    ])


    # Parse BLAST XML results
    with open(xml_file) as f:
        records = NCBIXML.parse(f)
        blast_record = next(records)


    # Build amino acid counts at each position
    seq_length = len(seq)

    position_counts = [{} for _ in range(seq_length)]


    # Count the query letter itself at each position
    for i, letter in enumerate(seq):
        position_counts[i][letter] = 1


    # Go through each alignment
    for alignment in blast_record.alignments:

        for hsp in alignment.hsps:

            query_pos = hsp.query_start - 1  # Convert to 0-based indexing

            query_seq = hsp.query
            subject_seq = hsp.sbjct

            for i in range(len(query_seq)):

                if query_seq[i] == "-":
                    continue  # Skip gaps in query

                subject_aa = subject_seq[i]

                if subject_aa == "-":
                    query_pos += 1
                    continue  # Skip gaps in subject

                if query_pos < seq_length:

                    if subject_aa not in position_counts[query_pos]:
                        position_counts[query_pos][subject_aa] = 0

                    position_counts[query_pos][subject_aa] += 1

                query_pos += 1


    # Calculate E6 and E20 for each position
    for i in range(seq_length):

        e6 = calculate_e6(position_counts[i])
        e20 = calculate_e20(position_counts[i])

        all_results.append({
            "Protein": protein,
            "Chain": chain,
            "No.": i,
            "Res": seq[i],
            "E6": e6,
            "E20": e20
        })


    print(
        f"  Done - {len(blast_record.alignments)} alignments found"
    )

    # Save checkpoint every 10 proteins
    if (idx + 1) % 10 == 0:
        pd.DataFrame(all_results).to_csv(checkpoint_file, index=False)
        print(f"  Checkpoint saved ({len(completed) + idx + 1} proteins)")


# Save results
results_df = pd.DataFrame(all_results)

results_df.to_csv(
    "/Users/patrickmincher/Desktop/LUSTIG SPSS/"
    "Generate_E6_E20_VKBAT/entropy_rewrite/entropy_results.csv",
    index=False
)

print("All done!")
