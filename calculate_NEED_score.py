'''
Date: Nov-05-2025
Author: Chih-Yuan Chou

Input profiles:
    1. raw SGReP z-scores profile with labels (.txt)
        e.g., TCGA_pan-cancer_PCC0.5_sgrep.txt
    2. HiSBiM pathway z-score profile (.txt)
        e.g., a_set_of_706_core_524_non-EGs_knowledge-based_pathway_enrichment.txt
Output profiles:
    1. selected pathways (positive pathways: 2, negative pathways: 1, others: 0)
        e.g., a_set_of_706_core_524_non-EGs_knowledge-based_pathway_enrichment_selected_pathways.txt
    2. random selected genes in each negative gene set
        e.g., t1_RF_random_non-EGs.txt
    3. mean of selected pathway importance by random selected negative gene sets
        e.g., mean_of_pathway_importances_by_all_t_cutoff.txt
    4. re-quantified gene scores and z-scores
        e.g., observed_gene_scores_and_pvalues_quantified_by_selected_pathways.txt
'''
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import os
from random import shuffle
from scipy import stats
from sklearn.ensemble import RandomForestClassifier

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run SIEG model")
    parser.add_argument('-s', '--sgrep', type=str, required=True,
                        help='Input SGReP z-scores profile')
    parser.add_argument('-g', '--grep', type=str, required=True,
                        help='Input a group of genes\' pathway z-scores profile')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--t_start', type=int, default=1,
                        help='Start of t')
    parser.add_argument('--t_end', type=int, default=2,
                        help='End of t')
    parser.add_argument('--t_step', type=int, default=1,
                        help='Step of t')
    args = parser.parse_args()
    f_sgrep = args.sgrep
    f_hisbim = args.grep
    output_path = Path(args.output)
    t_start = args.t_start
    t_end = args.t_end
    t_step = args.t_step
    output_path.mkdir(parents=True, exist_ok=True)

def calculate_g_score(df_h,df_pi,para_name,output_colName):
    pos_list = []
    neg_list = []
    pos_list = df_h[df_h[para_name]==2]["Pathway"].to_list()
    neg_list = df_h[df_h[para_name]==1]["Pathway"].to_list()
    pos_score_sum = df_pi[df_pi.columns[df_pi.columns.isin(pos_list)]].sum(axis=1)
    neg_score_sum = df_pi[df_pi.columns[df_pi.columns.isin(neg_list)]].sum(axis=1)
    df = pd.DataFrame(pos_score_sum-neg_score_sum, columns=[output_colName])
    return df
def calculate_zscore(para,obs_s):
    mean = np.mean(para)
    std_dev = np.std(para)
    z = (obs_s-mean)/std_dev
    return z

df_sgrep = pd.read_csv(f_sgrep, sep="\t",index_col=0)
df_sgrep.dropna(how="all", axis=0, inplace=True)
df_sgrep_num = df_sgrep.drop(['label'], axis=1)
df_sgrep_num.columns = df_sgrep_num.columns.to_list()
df_hisbim = pd.read_csv(f_hisbim, sep="\t")
dfs_feat_imp, dfs_g_score, dfs_g_pvalue = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

# calculate t score
df_hisbim.columns = ["System", "Subsystem","Pathway","pos","neg"]
df_hisbim["t_score"] = df_hisbim["pos"]-df_hisbim["neg"]
for t in range(t_start,t_end,t_step):
    # select positive and negative pathways
    cond1 = (df_hisbim["t_score"]>=t)
    cond2 = (df_hisbim["t_score"]<=-t)
    colName_t = "t"+str(t)
    df_hisbim[colName_t] = np.select([cond1,cond2],[2,1], default=0)
    file_name_selected_pathways = os.path.basename(f_hisbim)[:-4]+"_selected_pathways.txt"
    df_hisbim.to_csv(output_path/file_name_selected_pathways, sep="\t", index=False)
    # calculate pathway importance by selected pathways for 50 times
    random_times = 50
    l_selected_pathways = df_hisbim[df_hisbim[colName_t]>0]["Pathway"].values
    df_RF_sgrep = pd.concat([df_sgrep["label"],df_sgrep[l_selected_pathways]], axis=1)
    df_pos = df_RF_sgrep[df_RF_sgrep["label"]==2]
    df_neg = df_RF_sgrep[df_RF_sgrep["label"]==1]
    num_pos_sample = df_pos.shape[0]
    df_random_neg_genes = pd.DataFrame(0,index=df_neg.index,columns=["random"+ str(c) for c in range(1,random_times+1)])
    l_df_training = []
    for r in range(random_times):
        df_random_neg = df_neg.sample(num_pos_sample)
        for i in df_random_neg.index:
            df_random_neg_genes.loc[i,"random"+str(r+1)] = 1
        df = pd.concat([df_pos,df_random_neg])
        l_df_training.append(df)
    rf_model = RandomForestClassifier(
        bootstrap = False , 
        class_weight= None ,
        max_depth = None, 
        max_features = 'sqrt', 
        min_samples_leaf = 1, 
        min_samples_split = 2,
        n_estimators= 500,
        random_state = 42 )
    df_feat_imp = pd.DataFrame()
    times = 1
    for df_training in l_df_training:
        X_train = df_training.drop(['label'], axis=1)
        y_train = df_training.label
        features_id = list(X_train.columns)
        rf_model.fit(X_train, y_train)
        feat_imp = pd.DataFrame({"random"+str(times): rf_model.feature_importances_}, index=features_id)
        df_feat_imp = feat_imp if df_feat_imp.empty else df_feat_imp.merge(feat_imp, left_index=True, right_index=True)
        times += 1
    dfs_feat_imp[colName_t] = df_feat_imp.mean(axis=1)
    dic_mean_imp = df_feat_imp.mean(axis=1).to_dict()
    # calculate final pathway regulation*importance score
    df_pathImp = df_sgrep[l_selected_pathways].apply(lambda x: x*dic_mean_imp[x.name])
    file_name_weight = colName_t+"_selected_pathway_weighted_scores.txt"
    df_pathImp.to_csv(output_path/file_name_weight, sep="\t", index=True)
    # quantify all gene scores
    df_gscore = calculate_g_score(df_hisbim,df_pathImp,colName_t,colName_t)
    dfs_g_score = df_gscore if dfs_g_score.empty else dfs_g_score.merge(df_gscore, left_index=True, right_index=True)
    # transfer gene scores into z-scores
    shuffle(l_selected_pathways)
    df_shuffle_pathImp = df_pathImp
    df_shuffle_pathImp.columns = l_selected_pathways
    df_shuffle_gscore = calculate_g_score(df_hisbim,df_shuffle_pathImp,colName_t,colName_t)
    df_pvalue = df_gscore[colName_t].apply(lambda x: calculate_zscore(df_shuffle_gscore[colName_t],x))
    df_pvalue.columns = [colName_t]
    max_value = df_pvalue.replace([np.inf,-np.inf],np.nan).max()
    df_pvalue.replace([np.inf,-np.inf],max_value,inplace=True)
    dfs_g_pvalue = pd.DataFrame(df_pvalue) if dfs_g_pvalue.empty else dfs_g_pvalue.merge(pd.DataFrame(df_pvalue), left_index=True, right_index=True)
file_name_random = colName_t+"_RF_random_notCEGs.txt"
df_random_neg_genes.to_csv(output_path/file_name_random, sep="\t")
dfs_feat_imp.to_csv(output_path/"mean_of_pathway_importances_by_all_t_cutoff.txt", sep="\t", index=True)
dfs_g_score.columns = [c + "_gene_score" for c in dfs_g_score.columns]
dfs_g_pvalue.columns = [c + "_SIEG_score" for c in dfs_g_pvalue.columns]
dfs_gene_all_scores = dfs_g_score.merge(dfs_g_pvalue, left_index=True, right_index=True)
dfs_gene_all_scores.to_csv(output_path/"observed_gene_scores_and_pvalues_quantified_by_selected_pathways.txt", sep="\t", index=True)
