"""
Script for fine-tuning the RFC scored on precision, recall, and accuracy using a stratified K-fold cross-validation.
The cross-validation is stratified to keep the distributions of the train and test splits as similar as possible.
The best parameters are used to train the RFC on the entire dataset and is saved for use in ORCA.
"""

import gzip
import os
import pickle
import time
from typing import Callable

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (ConfusionMatrixDisplay, accuracy_score, classification_report, confusion_matrix,
                             make_scorer, precision_score, recall_score)
from sklearn.model_selection import (GridSearchCV, RandomizedSearchCV,
                                     StratifiedKFold, cross_validate,
                                     train_test_split)
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

RANDOM_STATE = 42
K = 5
PROCESSES = None  # os.cpu_count() - 1 # Number of processes started by the gridsearch

DATA_PATH_DORIC = "data/output/machine_learning/labels.csv"
DATA_PATH_EXP_SET = "data/output/machine_learning/labels_against_NCBI.csv"
MODEL_OUT_PATH = "data/output/machine_learning/ORCA_RFC_model_1_4_0.pkl"


def load_data_labels_from_path(path: str) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path)
    # !!! Keep this column order! This is the same order as ORCA will input the features when trying to predict.
    # ORCA uses a model trained on numpy array instead of a pandas dataframe.
    X = df[['Z_score', 'G_score', 'D_score']]

    # Feature-scaling does not provide any better performance
    # X_scaled = StandardScaler().fit_transform(X)
    y = df['correct']
    return X, y


def load_data_labels_from_df(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    X = df[['Z_score', 'G_score', 'D_score']]
    y = df['correct']
    return X, y


def validate_model(X: pd.DataFrame, y: pd.Series, classifier, k: int) -> dict:
    '''Cross-validates a given classifier in `k` folds based on accuracy, precision, and recall. Returns a dictionary with the results of each fold'''
    scoring = {
        'precision': make_scorer(precision_score),
        'recall': make_scorer(recall_score),
    }

    k_fold = StratifiedKFold(n_splits=k, shuffle=True, random_state=RANDOM_STATE)

    return cross_validate(estimator=classifier, X=X, y=y, cv=k_fold, scoring=scoring)


def fine_tune(X: pd.DataFrame, y: pd.Series, model, parameters: dict, n_jobs: int, refit_strategy: Callable, verbose: int = -1) -> GridSearchCV:
    cv = StratifiedKFold(n_splits=K, shuffle=True, random_state=RANDOM_STATE)
    start = time.perf_counter()
    grid_search = GridSearchCV(model, parameters, cv=cv, n_jobs=n_jobs, verbose=verbose,
                               scoring=["precision", "recall"], refit=refit_strategy).fit(X, y)
    # grid_search = RandomizedSearchCV(model, parameters, n_iter=10, cv=cv, n_jobs=n_jobs, verbose=verbose, scoring=["precision", "recall"], refit=refit_strategy, random_state=RANDOM_STATE).fit(X, y)
    print("TIME: ", time.perf_counter() - start)
    return grid_search


def print_CV_results(filtered_cv_results):
    """Pretty print for filtered CV results"""
    for mean_precision, std_precision, mean_recall, std_recall, params in zip(
        filtered_cv_results["mean_test_precision"],
        filtered_cv_results["std_test_precision"],
        filtered_cv_results["mean_test_recall"],
        filtered_cv_results["std_test_recall"],
        filtered_cv_results["params"],
    ):
        print(
            f"precision: {mean_precision:0.3f} (+/- {std_precision:0.03f}),"
            f" recall: {mean_recall:0.3f} (+/- {std_recall:0.03f}),"
            f" for {params}"
        )
    print()


def refit_strategy(cv_results: dict[str, np.ndarray]) -> int:
    """Define the strategy to select the best estimator.

    The strategy defined here is to filter-out all results below a precision threshold
    of 0.95, choose most performant model based on recall from these remaining ones.

    Parameters
    ----------
    cv_results : dict of numpy (masked) ndarrays
        CV results as returned by the `GridSearchCV`.

    Returns
    -------
    best_index : int
        The index of the best estimator as it appears in `cv_results`.
    """
    # print the info about the grid-search for the different scores
    precision_threshold = 0.8

    cv_results_ = pd.DataFrame(cv_results)
    print("All grid-search results:")
    print_CV_results(cv_results_)

    # Filter-out all results below the threshold
    high_precision_cv_results = cv_results_[
        cv_results_["mean_test_precision"] > precision_threshold
    ]

    print(f"Models with a precision higher than {precision_threshold}:")
    print_CV_results(high_precision_cv_results)

    high_precision_cv_results = high_precision_cv_results[
        [
            "mean_score_time",
            "mean_test_recall",
            "std_test_recall",
            "mean_test_precision",
            "std_test_precision",
            "rank_test_recall",
            "rank_test_precision",
            "params",
        ]
    ]

    # Select the most performant model from high_precision models in terms of recall
    print(
        "Out of the previously selected high precision models,\n"
        "we keep the model with the highest recall:\n\n"
        f"{high_precision_cv_results.loc[high_precision_cv_results['mean_test_recall'].idxmax()]}"
    )
    return high_precision_cv_results["mean_test_recall"].idxmax()


def main_tuning() -> None:
    X, y = load_data_labels_from_path(DATA_PATH_EXP_SET)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=RANDOM_STATE, stratify=y)

    model_parameters = {
        'bootstrap': [True, False],
        'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
        'max_features': ['log2', 'sqrt'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
    }

    # model_parameters = {
    #     'kernel': ['rbf', 'sigmoid'],
    #     'C': [10, 50, 100, 500, 1000, 1500, 2000],
    #     'gamma': ['scale', 'auto', 0.01, 0.05, 0.1, 0.15, 1]
    # }

    print("Start gridsearch")
    grid_search_model = fine_tune(
        X=X_train,
        y=y_train,
        model=RandomForestClassifier(random_state=RANDOM_STATE),  # SVC(random_state=RANDOM_STATE)
        parameters=model_parameters,
        refit_strategy=refit_strategy,
        n_jobs=PROCESSES,
        verbose=4
    )

    print()
    print("best parameters: ")
    for k, v in grid_search_model.best_params_.items():
        print(k, v)
    print()

    y_pred = grid_search_model.predict(X_test)
    print(classification_report(y_test, y_pred))


def main_DoriC_vs_exp_set() -> None:

    DoriC = pd.read_csv(DATA_PATH_DORIC)
    exp_set = pd.read_csv(DATA_PATH_EXP_SET)

    # only one version per accession in dataset, no need to check version
    DoriC_E = DoriC[~DoriC['accession'].isin(exp_set['accession'])]

    X_exp_set, y_exp_set = load_data_labels_from_df(exp_set)
    X_DoriC_E, y_DoriC_E = load_data_labels_from_df(DoriC_E)

    RFC_params = {
        'n_estimators': 600,
        'min_samples_split': 2,
        'min_samples_leaf': 2,
        'max_features': 'log2',
        'max_depth': 70,
        'bootstrap': False
    }

    print("Training on DoriC_E; testing on exp_set")
    RFC_model = RandomForestClassifier(random_state=RANDOM_STATE, **RFC_params)
    RFC_model.fit(X_DoriC_E, y_DoriC_E)
    y_pred_exp_set_by_DoriC_E = RFC_model.predict(X_exp_set)
    print(classification_report(y_exp_set, y_pred_exp_set_by_DoriC_E))
    cm = confusion_matrix(y_exp_set, y_pred_exp_set_by_DoriC_E, labels=RFC_model.classes_)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=RFC_model.classes_)
    disp.plot()
    plt.show()

    SVC_params = {
        'C': 50,
        'gamma': 0.1,
        'kernel': 'rbf'
    }

    print("Training on exp; testing on DoriC_E")
    SVC_model = SVC(random_state=RANDOM_STATE, **SVC_params)
    SVC_model.fit(X_exp_set, y_exp_set)
    y_pred_DoriC_E_by_exp_set = SVC_model.predict(X_DoriC_E)
    print(classification_report(y_DoriC_E, y_pred_DoriC_E_by_exp_set))
    cm = confusion_matrix(y_DoriC_E, y_pred_DoriC_E_by_exp_set, labels=SVC_model.classes_)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=SVC_model.classes_)
    disp.plot()
    plt.show()


def main_save_model() -> None:
    X, y = load_data_labels_from_path(DATA_PATH_DORIC)

    X = X.to_numpy()
    y = y.to_numpy()

    RFC_params = {
        'n_estimators': 600,
        'min_samples_split': 2,
        'min_samples_leaf': 2,
        'max_features': 'log2',
        'max_depth': 70,
        'bootstrap': False
    }

    RFC_model = RandomForestClassifier(random_state=RANDOM_STATE, **RFC_params)
    RFC_model.fit(X, y)
    with gzip.open(MODEL_OUT_PATH + ".gz", "wb") as fh:
        pickle.dump(RFC_model, fh)

    joblib.dump(RFC_model, MODEL_OUT_PATH)


if __name__ == '__main__':
    # main_tuning()
    # main_DoriC_vs_exp_set()
    main_save_model()
