"""This script is for measuring the performance of an SVC trained on dimensions provided by ORCA."""

from sklearn.svm import SVC

from sklearn.metrics import precision_score, recall_score, accuracy_score, make_scorer
from sklearn.model_selection import KFold, StratifiedKFold, GridSearchCV, cross_validate

import pandas as pd
import numpy as np

PATH = 'data\output\machine_learning\labels_against_NCBI.csv'
RANDOM_STATE = 42
K = 5
PROCESSES = 3 # Number of processes started by the gridsearch

def load_data(path):
    df = pd.read_csv(path)
    X = df[['Z_score', 'G_score', 'D_score', 'total_pot']]
    y = df['correct']
    return X, y


def validate_model(X, y, classifier, k):
    '''Cross-validates a given classifier in `k` folds based on accuracy, precision, and recall. Returns a dictionary with the results of each fold'''
    scoring = {
        'accuracy' : make_scorer(accuracy_score),
        'precision' : make_scorer(precision_score),
        'recall' : make_scorer(recall_score),
    }

    k_fold = StratifiedKFold(n_splits=k, shuffle=True, random_state=RANDOM_STATE)

    return cross_validate(estimator=classifier, X=X, y=y, cv=k_fold, scoring=scoring) 


def fine_tune(X, y, model, parameters, n_jobs, verbose=-1):
    cv = StratifiedKFold(n_splits=K, shuffle=True, random_state=RANDOM_STATE)
    grid_search = GridSearchCV(model, parameters, cv=cv, n_jobs=n_jobs, verbose=verbose, scoring=["precision", "recall", "accuracy"], refit="precision").fit(X, y)
    return grid_search.best_params_


def main():
    X, y = load_data(PATH)
    model_parameters = {
        'C': [0.001, 0.1, 1, 5, 10, 30, 50, 60, 80, 100],       # 'degree' parameter might be added later, depends on preffered kernel (only used by 'poly')
        # 'degree': [i for i in range(10)],
        'kernel': ['linear', 'poly', 'sigmoid', 'rbf']  # 'precomputed' only works on NxN datasets
    }

    best_params = fine_tune(X, y, SVC(random_state=RANDOM_STATE), model_parameters, PROCESSES, True)
    #                                                                       [params]
    SVC_scores_tuned   = validate_model(X, y, SVC(**best_params), K)
    SVC_scores_default = validate_model(X, y, SVC(), K)

    print(f'''
Fine-tuned :
precision
 mean : {np.mean(SVC_scores_tuned['test_precision']):.3f}
 std  : {np.std(SVC_scores_tuned['test_precision']):.3f} 
recall
  mean : {np.mean(SVC_scores_tuned['test_recall']):.3f}
  std  : {np.std(SVC_scores_tuned['test_recall']):.3f} 
accuracy
  mean : {np.mean(SVC_scores_tuned['test_accuracy']):.3f}
  std  : {np.std(SVC_scores_tuned['test_accuracy']):.3f}
    ''')
    print(f'''
Default :
precision
 mean : {np.mean(SVC_scores_default['test_precision']):.3f}
 std  : {np.std(SVC_scores_default['test_precision']):.3f} 
recall
  mean : {np.mean(SVC_scores_default['test_recall']):.3f}
  std  : {np.std(SVC_scores_default['test_recall']):.3f} 
accuracy
  mean : {np.mean(SVC_scores_default['test_accuracy']):.3f}
  std  : {np.std(SVC_scores_default['test_accuracy']):.3f}
    ''')

if __name__ == '__main__':
    main()