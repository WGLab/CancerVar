#coding=utf-8
# Zilin Ren
# 2021-01-17

import argparse, sys, re

import numpy  as np
import pandas as pd
from sklearn.impute import KNNImputer
import uuid 



def get_args():

    parser = argparse.ArgumentParser(description='feature creator from cancervar output')

    parser.add_argument( "-a", "--annovar_path", help = "the path to annovar file", \
        default=None, type=str, required=True )
    
    parser.add_argument( "-c", "--cancervar_path", help = "the path to cancervar file", \
        default=None, type=str, required=True )    

    parser.add_argument( "-m", "--method", help = "output evs features or ensemble features (option: evs, ensemble)", \
        default="ensemble", type=str, required=False )

    parser.add_argument( "-n", "--missing_count", help="variant with more than N missing features will be discarded, (default: 5)", \
        default=5, type=int, required=False )

    parser.add_argument( "-d", "--database", help="database for feature normalization", \
        default = None, type=str, required=True)

    parser.add_argument( "-o", "--output", help="the path to output", \
        default=None, type=str, required=True )

    ## run
    args = parser.parse_args()
    return args



def main():
    args = get_args()
    
    ## anno feature names
    fun_name = get_continuous_feature()

    ## evs feature names
    evs_name = get_evs_features()

    ## gaussian params
    gaussian_params = [False, 0, 0.02]

    evs_feature, fun_feature = data_loader(args, fun_name, evs_name, gaussian_params)

    if args.method == "evs":
        data_saver(evs_feature, args.output)

    else:
        filter_feature_dat = data_filter_imputer(evs_feature, fun_feature, args.missing_count, args.database)
        data_saver(filter_feature_dat, args.output)

    return 1


##
def data_saver(dat, output):
    dat.to_csv(output, header = None, index=True,sep="\t")
    return 1


##
def data_filter_imputer(evs_feature, fun_feature, missn, database_path):
    filtered_fun = fun_feature.loc[fun_feature.isna().sum(axis=1) <= missn, :]

    if filtered_fun.shape[0] == 0:
        return filtered_fun

    filtered_evs = evs_feature.loc[filtered_fun.index, :]

    db_dat      = np.load(database_path, allow_pickle=True)
    tot_dat     = np.vstack( (db_dat, filtered_fun) )
    for col in range(tot_dat.shape[1]):
        tot_dat[:, col] = ( tot_dat[:, col] - min(tot_dat[:, col]) ) / ( max(tot_dat[:, col]) - min(tot_dat[:, col]) )

    if missn > 0:
        imputer     = KNNImputer(n_neighbors=40)
        tot_dat = imputer.fit_transform(tot_dat)


    filled_feature = tot_dat[db_dat.shape[0]:]
    filled_df      = pd.DataFrame(filled_feature, columns = filtered_fun.columns)
    filled_df.set_index(filtered_evs.index, inplace=True)
    
    merged_dat = pd.concat([filled_df, filtered_evs], axis=1, ignore_index=True)

    return merged_dat



def data_loader(args, fun_name, evs_name, gaussian_params):

    cancervar_path, annovar_path = args.cancervar_path, args.annovar_path

    cancervar_dat = pd.read_csv(cancervar_path, sep='\t', dtype='str', na_values=".")
    cancervar_dat['index'] = cancervar_dat.apply(lambda x: 'row_{%s}_chr{%s}_start{%s}_end{%s}_ref{%s}_alt{%s}' % (x.name, x['#Chr'], x['Start'], x["End"], x['Ref'], x['Alt']), axis = 1)

    annovar_dat = pd.read_csv(annovar_path, sep='\t', na_values=".")    
    annovar_dat['index']   = annovar_dat.apply(lambda x: 'row_{%s}_chr{%s}_start{%s}_end{%s}_ref{%s}_alt{%s}' % (x.name, x['Chr'], x['Start'], x["End"], x['Ref'], x['Alt']), axis = 1)

    evs_feature = func_individual_evs_dat(cancervar_dat, evs_name, gaussian_params)

    if args.method == "evs":
        return evs_feature, None
    
    if not all( cancervar_dat['index'] == annovar_dat['index'] ):
        print("ERROR! the row of %s is different with %s " % (annovar_path, cancervar_path))
        sys.exit(1)
    
    fun_feature = func_individual_function_dat(cancervar_dat, annovar_dat, fun_name)
    
    return evs_feature, fun_feature


## 
def get_continuous_feature():
    features_with_continuous_val = ["SIFT_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "LRT_score", "MutationTaster_score", 
                                    "MutationAssessor_score", "FATHMM_score", "PROVEAN_score", "VEST3_score", "CADD_raw", 
                                    "CADD_phred", "DANN_score", "fathmm-MKL_coding_score", "MetaSVM_score", "MetaLR_score", 
                                    "integrated_fitCons_score","integrated_confidence_value","GERP++_RS", "phyloP7way_vertebrate", "phyloP20way_mammalian", 
                                    "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", ]
    return features_with_continuous_val


## 
def get_evs_features():
    evs_val  = ["EVS_1", "EVS_2", "EVS_3", "EVS_4", "EVS_5", "EVS_6", "EVS_7", "EVS_8", "EVS_9", "EVS_10", "EVS_11", "EVS_12"]
    colnames = ["%s_{%s}" % (i, j) for i in evs_val for j in [-1, 0, 1, 2]]
    return colnames



## 
def func_individual_function_dat(cancervar_dat, annovar_dat, fun_name):
    feature_list = []

    for colname in fun_name:
        # colname_list.append(colname)
        if colname in  cancervar_dat.columns:
            feature_list.append( cancervar_dat.loc[:, [colname]] )
        else:
            feature_list.append( annovar_dat.loc[:, [colname]] )

    function_features = pd.concat(feature_list, axis=1)
    Index    = cancervar_dat['index']
    return function_features.set_index(Index).astype('float')


## 
def func_individual_evs_dat(cancervar_dat, evs_name, params):
    # params = (False, 0, 0.02)
    def get_evidence_feture(string, params):

        add_gaussian, mu, sigma = params

        evidence_string  = re.search(r'=(\[.*\])', string).group(1)
        evidence_feature = eval(evidence_string)

        dummy_feature = []
        if add_gaussian:
            gaussian_noise = np.random.normal(mu, sigma, 4)
            for ind, indice in enumerate(evidence_feature):
                if indice == -1:
                    dummy_feature += list(np.array([1, 0, 0, 0]) + gaussian_noise)
                elif indice == 0:
                    dummy_feature += list(np.array([0, 1, 0, 0]) + gaussian_noise)
                elif indice == 1:
                    dummy_feature += list(np.array([0, 0, 1, 0]) + gaussian_noise)
                else:
                    dummy_feature += list(np.array([0, 0, 0, 1]) + gaussian_noise)
        else:
            for ind, indice in enumerate(evidence_feature):
                # -1, 0, 1, 2
                if indice == -1:
                    dummy_feature += [1, 0, 0, 0] 
                elif indice == 0:
                    dummy_feature += [0, 1, 0, 0]
                elif indice == 1:
                    dummy_feature += [0, 0, 1, 0]
                else:
                    dummy_feature += [0, 0, 0, 1]

        return  dummy_feature


    Index    = cancervar_dat['index']
    
    evidence_series  = cancervar_dat[" CancerVar: CancerVar and Evidence "].apply(lambda x: get_evidence_feture(x, params))
    evidence_feature = pd.DataFrame(evidence_series.tolist(), columns=evs_name)
    return evidence_feature.set_index(Index)


if __name__ == '__main__':
    main()
