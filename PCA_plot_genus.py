import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re

file_path = '/workdir/ls764/PCA/all_samples'
pc_path ="/workdir/ls764/PCA/all_samples"

def sort_by_state(surface):
    """
    Sort the samples by the surface type 
    Return a nested list which the first element is a list that contains all the sample dataframe and the second element contain a 
    list that stores the sample name and the third element contains a list that stores its corresponding state
    """
    sample_name_table = pd.read_csv('sample_name.csv',header=None)
    sample_name_table.columns = ["sample_name","state"]
    motu_table_path = '/workdir/ls764/PCA/motus_outputs/{sample_name}/{sample_name}.motus.ratio'
    motus = []
    #sample_name_table = sample_name_table[sample_name_table["state"].str.contains(surface)] #exlude control samples

    sample_name_table = sample_name_table[sample_name_table["state"].str.contains('Control')|sample_name_table["state"].str.contains(surface)]
    #Exclude outlier
    # if surface == "implant":
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("I_11_MB")]
    #     #sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("IV_1_MB")] #2x
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("IV_1_ML")]
    #     sample_name_table = sample_name_table[~sample_name_table['state'].str.contains("blank")]
    # if surface == "tooth":
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("2_31_MB")]
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("1_31_ML")]#2x
    #     sample_name_table = sample_name_table[~sample_name_table['state'].str.contains("blank")]
    # if surface == " ":
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("2_31_MB")]
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("I_11_MB")]
    #     sample_name_table = sample_name_table[~sample_name_table['state'].str.contains("blank")]
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("IV_1_ML")]
    #     sample_name_table = sample_name_table[~sample_name_table['sample_name'].str.contains("IV_1_MB")] #2x

    for sn in sample_name_table['sample_name'].values:
        motu_table_fn = motu_table_path.format(sample_name=sn)
        motus.append(pd.read_csv(motu_table_fn, sep='\t',header=None))
    samples =[]
    samples.append(motus)
    sample_name = sample_name_table['sample_name'].tolist()
    state = sample_name_table['state'].tolist()
    samples.append(sample_name)
    samples.append(state)
    return samples

def get_genus(df,sample_name):
    """
    Filter out the genus abundance from the filepath and create a dataframe with a column name of sample_name_Abundance
    """
    #df = pd.read_csv(filepath,delimiter='\t',header=None)
    df = df.drop(0)
    df.columns = ["Taxon","Abundance"]
    total_ab = float(df.loc[1,"Abundance"])
    pipe_count = df['Taxon'].str.count('\|')
    filtered_genus = df[pipe_count == 5] #filtered all genus level
    filtered_genus["Abundance"] = filtered_genus["Abundance"].astype(float)
    filtered_genus['Taxon'] = filtered_genus['Taxon'].str.rsplit('|', 1).str[-1] #extract the name of the genus 
    incertae_sedis = filtered_genus[filtered_genus['Taxon'].str.contains("incertae sedis")] #extract the unknow species percentage
    uncertain_ab = list(incertae_sedis["Abundance"])
    drop_index = incertae_sedis.index.tolist()
    filtered_genus.drop(drop_index,axis=0,inplace= True)
    filtered_genus[sample_name+'Relative_ab'] = filtered_genus["Abundance"].div(total_ab)
    filtered_genus = filtered_genus.groupby('Taxon')[sample_name+'Relative_ab'].sum().reset_index()
    print(filtered_genus[sample_name+'Relative_ab'].sum())
    print
    return filtered_genus

def merge(samples):
    """
    Merge all the sample dataframe into a formatted dataframe that can be used for PCA
    Take in a list that contain all the samples info
    """
    #samples_name = ["Sample" + str(i+1) for i in range(len(samples))]
    null = pd.DataFrame({'Taxon': []})
    for i in range(len(samples[0])):
        # filepath = samples[i]
        name = samples[1][i]       
        current = get_genus(samples[0][i],name)
        null = pd.merge(null, current, on='Taxon', how='outer')   
    df = null.T
    df = df.fillna(0)
    return df

def pca_plot(df,samples,surface):
    Taxon = df.loc['Taxon']
    df.columns = Taxon
    df = df.drop(labels='Taxon')
    df['Labels']= samples[2]
    df.to_csv('{}/merged_samples{}.csv'.format(file_path,surface))
    x = df.loc[:, Taxon].values
    x = StandardScaler().fit_transform(x) # normalizing the abundance
    # Taxon_cols = ['Taxon'+str(i) for i in range(x.shape[1])]
    # normalised_taxon = pd.DataFrame(x,columns=Taxon_cols)
    #create pca
    pca_taxon = PCA(n_components=2)
    # pca_taxon = PCA()
    principalComponents_taxon = pca_taxon.fit_transform(x)
    per_var = np.round(pca_taxon.explained_variance_ratio_*100,decimals=1) #Importance of each PC    
    print('Explained variation per principal component: {}'.format(pca_taxon.explained_variance_ratio_))
    principal_taxon_Df = pd.DataFrame(data = principalComponents_taxon, columns = ['principal component 1', 'principal component 2'])
    principal_taxon_Df['sample name'] = samples[1]
    principal_taxon_Df['state'] = samples[2]
    # principal_taxon_Df.to_csv('{}/PC_values_{}'.format(file_path,surface))
    plt.figure()
    plt.figure(figsize=(10,10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    plt.xlabel('PC1 - {0}%'.format(per_var[0]),fontsize=20)
    plt.ylabel('PC2 - {0}%'.format(per_var[1]),fontsize=20)
    plt.title("Principal Component Analysis of Bacteria distribution_"+surface,fontsize=20)
    targets = np.unique(samples[2])
    colors=['r', 'g','b','c','m','y','k','tab:purple','tab:pink']
    # colors = []
    # for i in range(len(colors_list)):
    #     colors.append(colors_list[i])
    df= df.reset_index()
    for target, color in zip(targets,colors):
        indicesToKeep = df['Labels'] == target
        plt.scatter(principal_taxon_Df.loc[indicesToKeep, 'principal component 1']
                , principal_taxon_Df.loc[indicesToKeep, 'principal component 2'], c = color, s = 50)
    plt.legend(targets,prop={'size': 15})
    plt.savefig("{}/Principal Component Analysis of Bacteria distribution_{}.png".format(file_path,surface))
    return x

def scree_plot(normalized_data,surface,samples,df):
    Taxon = df.loc['Taxon']
    plt.figure()
    pca_taxon = PCA()
    principalComponents_taxon = pca_taxon.fit_transform(normalized_data)
    for i in range(len(pca_taxon.components_)):
        loading_scores = pd.Series(pca_taxon.components_[i],index=Taxon)
        sorted_loading_scores= loading_scores.abs().sort_values(ascending=False)
        top_15_genus = sorted_loading_scores[0:15].index.values
        DF = sorted_loading_scores.reset_index()
        DF.to_csv("{}/top_10_genus_PCA{}.csv_{}".format(pc_path,i+1,surface),index=False, header=False)
    cols = ["principal component " + str(i+1) for i in range(len(principalComponents_taxon))]
    principal_taxon_Df = pd.DataFrame(data = principalComponents_taxon, columns = cols)
    principal_taxon_Df['sample name'] = samples[1]
    principal_taxon_Df['state'] = samples[2]
    principal_taxon_Df.to_csv('{}/PC_values_{}'.format(file_path,surface))
    per_var = np.round(pca_taxon.explained_variance_ratio_*100,decimals=1) 
    labels = ['PC'+str(x) for x in range (1,len(per_var)+1)]
    plt.bar(x=range(1,len(per_var)+1),height=per_var,tick_label= labels)
    plt.xlabel('PC')
    plt.ylabel("Percentage of Explained Variance")
    plt.title('Scree plot')
    plt.savefig("{}/Principal Component Scree plot_{}.png".format(file_path,surface))


def main(): 
    # all_sam = sort_by_state(' ')
    # # all_df =merge(all_sam)
    # # x_all=pca_plot(all_df,all_sam,"all")
    implant_sam = sort_by_state('implant')
    implant_df =merge(implant_sam)
    x_implant=pca_plot(implant_df,implant_sam,"implant")
    # tooth_sam = sort_by_state('tooth')
    # tooth_df =merge(tooth_sam)
    # x_tooth = pca_plot(tooth_df,tooth_sam,"tooth")
    # #scree_plot(x_all,"all",all_sam,all_df)
    scree_plot(x_implant,"implant",implant_sam,implant_df)
    # scree_plot(x_tooth,"tooth",tooth_sam,tooth_df)
    # df = get_genus(all_sam[0][1],"1")
    # print (df)




    # together = sample_dict
    # tooth = {"Stage 3","Healthy T","Stage 1",'Blank',"Control E. coli and L. plantarum"}
    # tooth_dict = {key: value for key, value in sample_dict.items() if value in tooth}
    # implant = {"Healthy I","Moderate","Mucositis","Severe","Control E. coli and L. plantarum","Blank"}
    # implant_dict={key: value for key, value in sample_dict.items() if value in implant}
    # pca_plot(sample_dict,"all")
    #pca_plot(tooth_dict,"tooth")
    #pca_plot(implant_dict,"implant")
    
if __name__ == '__main__':
    main()








