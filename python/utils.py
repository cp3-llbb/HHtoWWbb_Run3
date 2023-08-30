import pandas as pd

def checkBatches(df,label_column,batch_size=128,weight_column='weight'):
    N_checks = 200
    labels = pd.unique(df[label_column])
    sum_label = {label:0. for label in labels}    
    N_label = {label:0 for label in labels}    
    for i in range(N_checks):
        rnd_df = df.sample(batch_size)
        for label in labels:
            sum_label[label] += rnd_df[rnd_df[label_column]==label][weight_column].sum()
            N_label[label] += rnd_df[rnd_df[label_column]==label][weight_column].shape[0]

    print (f'On average, per batch the total weight is')
    for label in labels:
        sum_label[label] /= N_checks
        N_label[label] /= N_checks
        print (f'\t... {label:20s}: {sum_label[label]:15.9f} [{N_label[label]} events]')
