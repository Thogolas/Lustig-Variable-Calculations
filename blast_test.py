import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML

df = pd.read_csv("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/entropy_rewrite/sequences.csv")
first_seq = df["Sequence"].iloc[0]

print(first_seq)

result_handle = NCBIWWW.qblast("blastp", "nr", first_seq)

with open("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/entropy_rewrite/test_blast.xml", "w") as f:
    f.write(result_handle.read())
print("Done")
