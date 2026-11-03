import pandas as pd


#Reading bondugula
df = pd.read_csv("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/merged.csv")

grouped = df.groupby(['Protein', 'Chain'])

#Each group takes the Res column and adds all letters into a single string
sequences = grouped['Res'].apply(''.join)

sequences = sequences.reset_index()
sequences.columns = ['Protein', 'Chain', 'Sequence']
sequences.to_csv("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/entropy_rewrite/sequences.csv", index=False)


