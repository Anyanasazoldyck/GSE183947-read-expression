#This script is to manipulate gene expression data
#"C:\Users\lenovo\Desktop\crpo\journal club"
# Import libraries
import pandas as pd
import numpy as np
import GEOparse





#read the data 
dat=pd.read_csv("./GSE183947_fpkm.csv")


# Get metadata

gse = GEOparse.get_GEO('GSE183947') #using GEOparse.get_GEO(GSE183947)will download a soft file into your dic
metadata = gse.metadata
metadata_dict=gse.phenotype_data #Get the phenotype data for each of the sample.

#metadata to dataframe
metadata_df = pd.DataFrame(metadata_dict)
metadata_df['tissue'] = metadata_df['characteristics_ch1.0.tissue'].str.replace("tissue: ", "") #for easier readability 
metadata_df['metastasis'] = metadata_df['characteristics_ch1.1.metastasis'].str.replace("metastasis: ", "")
metadata_subset =  metadata_df[['title', 'tissue', 'metastasis', 'description']]
metadata_subset.to_csv("./metadata_clean.csv")




#reshape the gse data 
dat.rename(columns={'Unnamed: 0': "gene"}, inplace=True)
dat_long = dat.melt(id_vars="gene", var_name="samples", value_name="FPKM")
print(dat_long.head(5))

# Join dataframes
dat_long = dat_long.merge(metadata_df, left_on="samples", right_on="description", how="left")
print(dat_long.head(5))

# Filter, group_by, summarize and arrange
result = (
    dat_long[dat_long['gene'].isin(['BRCA1', 'BRCA2'])]
    .groupby(['gene', 'tissue'])
    .agg(mean_FPKM=('FPKM', 'mean'), median_FPKM=('FPKM', 'median'))
    .sort_values('mean_FPKM', ascending=False)
)

print(result)



