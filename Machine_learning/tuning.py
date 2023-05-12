from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression

from sklearn.metrics import precision_score
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV

import pandas as pd
import numpy as np

PATH = 'Machine_Learning/tuning.csv'
RANDOM_STATE = 42
K = 5
PROCESSES = 3 # Number of processes started by the gridsearch

def load_data(path):
    df = pd.read_csv(path)
    X = df[['Z_occurance', 'G_occurance', 'D_occurance']]
    y = df['Correct']
    return X, y


def train_model(X, y, classifier, k):
    '''input the classifier, the training data and corresponding labels and returns a trained classifier and the accuracy of a prediction with k-fold validation'''
    
    # initialise the model and k-fold cross validator
    model = classifier
    k_fold = KFold(n_splits=k)
    acc_score = []

    # split the data into k folds
    for train_index, vali_index in k_fold.split(X):
        
        # assign the training and validation data accordingly
        train_data, vali_data = X[train_index, :], X[vali_index, :]
        train_labels, vali_labels = y[train_index], y[vali_index]

        # fit the model onto the training data
        model.fit(train_data, train_labels)

        # predict the validation data
        prediction = model.predict(vali_data)

        # assess the accuracy of the prediction
        accuracy_of_current_fold = precision_score(prediction , vali_labels)
        acc_score.append(accuracy_of_current_fold)

    # average all the accuracies
    avg_acc_score = sum(acc_score)/k

    return model, avg_acc_score


def fine_tune(X, y, model_dict, parameter_dict, n_jobs, print_found_parameters=True):
    best_params_models = []
    for name, parameters in parameter_dict.items():
        model = model_dict[name]
        cv = KFold(n_splits=K)
        grid_search = GridSearchCV(model, parameters, cv=cv, n_jobs=n_jobs, verbose=4, scoring="precision").fit(X, y)
        best_params_models.append([name, grid_search.best_params_])

    if print_found_parameters:
        for i in best_params_models:
            print(i)
    return best_params_models


if __name__ == '__main__':
    X, y = load_data(PATH)
    model_parameters = {
        "SVC": {
            'C': [1, 10, 100, 1000],                        # 'degree' parameter might be added later, depends on preffered kernel (only used by 'poly')
            'kernel': ['linear', 'poly', 'sigmoid', 'rbf']  # 'precomputed' only works on NxN datasets
        },
        "LogisticRegression": {
            'penalty': ['l1', 'l2', 'elasticnet', 'none'],
            'C': [1, 10, 100, 1000],
            'max_iter': [x for x in range(100, 600, 100)],
            'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga']
        }
    }

    models = {
        "SVC": SVC(random_state=RANDOM_STATE),
        "LogisticRegression": LogisticRegression(random_state=RANDOM_STATE)
    }

    best_params = fine_tune(X, y, models, model_parameters, PROCESSES)
    #                                                             [name][params]
    best_SVC_model_best, SVC_acc_best = train_model(SVC(**best_params[0][1]), X, y, K)
    best_SVC_model_def, SVC_acc_def = train_model(SVC(), X, y, K)
    best_LoR_model_best, LoR_acc_best = train_model(LogisticRegression(**best_params[1][1]), X, y, K)
    best_LoR_model_def, LoR_acc_def = train_model(LogisticRegression(), X, y, K)

# [CV 3/3] END .................C=100, kernel=rbf;, score=0.944 total time=  21.8s

    print(f'SVC_finetune    : {SVC_acc_best:.3f}')
    print(f'SVC_default     : {SVC_acc_def:.3f}')
    print(f'LogReg_finetune : {LoR_acc_best:.3f}')
    print(f'LogReg_default  : {LoR_acc_def:.3f}')
