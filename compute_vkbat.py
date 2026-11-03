import pandas as pd
from entropy_calc import calculate_vkbat

df = pd.read_csv("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/merged.csv")

predictor_cols = ['chou_fasman', 'sspro_5', 'gor4', 'dsc', 'jnet', 'psipred']
df['new_vkbat'] = df[predictor_cols].apply(lambda row:calculate_vkbat(row.tolist()), axis=1)

df.to_csv("/Users/patrickmincher/Desktop/LUSTIG SPSS/Generate_E6_E20_VKBAT/entropy_rewrite/merged_with_vkbat.csv", index = False)
print("Done! Rows:", len(df))