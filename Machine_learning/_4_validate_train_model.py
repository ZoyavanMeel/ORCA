"""
Script for fine-tuning the SVM scored on precision, recall, and accuracy using a stratified K-fold cross-validation.
The cross-validation is stratified to keep the distributions of the train and test splits as similar as possible.
The best parameters are used to train the SVC on the entire dataset and is saved for use in ORCA.
"""

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (accuracy_score, make_scorer, precision_score,
                             recall_score)
from sklearn.model_selection import (GridSearchCV, StratifiedKFold,
                                     cross_validate)
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

RANDOM_STATE = 42
K            = 5
PROCESSES    = 32 # Number of processes started by the gridsearch

CLUSTER        = "/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/ORCA/"
DATA_PATH      = CLUSTER + "data/output/machine_learning/labels.csv"
MODEL_OUT_PATH = CLUSTER + "data/output/machine_learning/24k_set_model.pkl"


def load_data_labels(path: str) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path)
    X = df[['Z_score', 'G_score', 'D_score', 'total_pot']]
    
    X_scaled = StandardScaler().fit_transform(X)
    y = df['correct']
    return X_scaled, y


def validate_model(X: pd.DataFrame, y: pd.Series, classifier, k: int) -> dict:
    '''Cross-validates a given classifier in `k` folds based on accuracy, precision, and recall. Returns a dictionary with the results of each fold'''
    scoring = {
        'accuracy' : make_scorer(accuracy_score),
        'precision' : make_scorer(precision_score),
        'recall' : make_scorer(recall_score),
    }

    k_fold = StratifiedKFold(n_splits=k, shuffle=True, random_state=RANDOM_STATE)

    return cross_validate(estimator=classifier, X=X, y=y, cv=k_fold, scoring=scoring) 


def fine_tune(X: pd.DataFrame, y: pd.Series, model, parameters: dict, n_jobs: int, verbose: int = -1) -> dict:
    cv = StratifiedKFold(n_splits=K, shuffle=True, random_state=RANDOM_STATE)
    grid_search = GridSearchCV(model, parameters, cv=cv, n_jobs=n_jobs, verbose=verbose, scoring=["precision", "recall", "accuracy"], refit="precision").fit(X, y)
    return grid_search.best_params_


def main() -> None:
    X, y = load_data_labels(DATA_PATH)
    model_parameters = {
        'C': [0.001, 0.1, 1, 5, 10, 30, 50, 60, 80, 100, 500, 1000], # 'degree' parameter might be added later, depends on preffered kernel (only used by 'poly')
        # 'degree': [i for i in range(10)],
        'gamma': ['scale', 'auto', 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000],
        'kernel': ['linear', 'poly', 'sigmoid', 'rbf']  # 'precomputed' only works on NxN datasets
    }

    print("Start gridsearch")
    best_params = fine_tune(
        X=X,
        y=y,
        model=SVC(random_state=RANDOM_STATE),
        parameters=model_parameters,
        n_jobs=PROCESSES,
        verbose=1
    )

    print()
    print("best parameters: ")
    for k, v in best_params.items():
        print(k, v)
    print()

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


    print("Training on full dataset...", end="")
    model = SVC(**best_params).fit(X, y)
    joblib.dump(model, MODEL_OUT_PATH)
    print("Done!")

if __name__ == '__main__':
    main()
