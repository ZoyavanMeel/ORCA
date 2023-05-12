import pandas as pd
import numpy as np
from sklearn.svm import SVC
import joblib

exp_refseq = [ # Accessions that have been experimentally verified.
    'NC_000964', 'NC_002947', 'NC_003272', 'NC_003869', 'NC_005090', 'NC_006461', 'NC_007633', 'NC_000913', 'NC_003047',
    'NC_007604', 'NC_000962', 'NC_002696', 'NC_002971', 'NC_005363', 'NC_008255', 'NC_009850', 'NC_010546', 'NC_010547', 'NC_011916'
]

all_samples_df = pd.read_csv('Hyperparameter tuning/tuning.csv')
raw_refseqs_df = pd.read_csv('DoriC data prep/DoriC_oriC_concat_entries.csv')['RefSeq']

# train on all but the experimental set
train_refseqs = raw_refseqs_df[~raw_refseqs_df.isin(exp_refseq)]
test_refseqs = raw_refseqs_df[raw_refseqs_df.isin(exp_refseq)]

train_samples_refseqs = []
test_samples_refseqs = []
for i, sample in all_samples_df.iterrows():
    if sample['RefSeq_oriC'][:-2] not in test_refseqs.to_list():
        train_samples_refseqs.append(sample['RefSeq_oriC'])
    else:
        test_samples_refseqs.append(sample['RefSeq_oriC'])

# .to_numpy() to get rid of feature name warning
X_train = all_samples_df[all_samples_df['RefSeq_oriC'].isin(train_samples_refseqs)][['Z_occurance', 'G_occurance', 'D_occurance']].to_numpy()
y_train = all_samples_df[all_samples_df['RefSeq_oriC'].isin(train_samples_refseqs)]['Correct'].to_numpy()

X_test = all_samples_df[all_samples_df['RefSeq_oriC'].isin(test_samples_refseqs)][['Z_occurance', 'G_occurance', 'D_occurance']].to_numpy()
y_test = all_samples_df[all_samples_df['RefSeq_oriC'].isin(test_samples_refseqs)]['Correct'].to_numpy()

# parameters from tuning.py
model = SVC(C=500, kernel='rbf', random_state=42).fit(X_train, y_train)

# save model as pickle-file
joblib.dump(model, 'exp_train_model.pkl')
