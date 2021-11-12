import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold
from interpret.glassbox import ExplainableBoostingClassifier
from pathlib import Path
import pickle
import os


def train_model_ebm(df_ensemble,
                    features, features_train=None, features_test=None,
                    chosen_model=ExplainableBoostingClassifier(),
                    feature_to_filter_na='delta_psi',
                    index_GroupKFold='sample', nsplits=5,
                    save_dir=None, write_to_pickle=False, save_results=False):

    results_all = list()
    models_all = list()

    if features_train is None:
        features_train = features
        features_test = features

    if feature_to_filter_na is not None:
        df = df_ensemble[
            ~df_ensemble[feature_to_filter_na].isna()
            ][[*features, 'outlier']]
        df_missing = df_ensemble[
            df_ensemble[feature_to_filter_na].isna()
            ][[*features, 'outlier']]
        X_missing = df_missing[features].fillna(0)
        y_missing = df_missing[['outlier']]
        y_pred_missing_iterm = np.zeros([df_missing.shape[0],nsplits])
    else:
        df = df_ensemble
        df_missing = None

    X = df[features]
    y = df[['outlier']]
    
    groups = df.index.get_level_values(index_GroupKFold)
    gkf = GroupKFold(n_splits=nsplits)
    for fold, (train, test) in enumerate(gkf.split(X, y, groups=groups)):
        # train, test
        X_train = X[features_train].iloc[train].fillna(0)
        X_test = X[features_test].iloc[test].fillna(0)
        y_train = y.iloc[train]
        y_test = y.iloc[test]
        # fit model
        model = ExplainableBoostingClassifier() #TODO: model=chosen_model
        model.fit(X_train, y_train)
        models_all.append(model)
        y_pred = model.predict_proba(X_test)[:, 1]  
        if write_to_pickle == True:
            _write_to_pickle(model, save_dir, fold)
        # store predictions of fold
        results = _update_fold_results(X, y, y_test, y_pred, test, fold, \
            features, features_train, features_test, model)    
        results_all.append(results)
    results_all_df = pd.concat(results_all) 
    
    # Join with outliers that we do not have predictions for (based on feature_to_filter_na)
    if df_missing is not None: 
        results_all_df = _update_results_with_missing(X_missing, y_missing, y_pred_missing_iterm,\
            features, features_train, features_test, models_all, results_all_df)
    
    if save_results == True:
        _save_results(results_all_df, save_dir)
    
    return results_all_df, models_all


def _write_to_pickle(model, save_dir, fold):
    path = Path(save_dir)
    path.mkdir(parents=True, exist_ok=True)
    pickle_filename = os.path.join(save_dir, 'cross_val=' + str(fold) + '.pkl')
    with open(pickle_filename, 'wb') as file:
        pickle.dump(model, file)


def _save_results(results_all_df, save_dir):
    path = Path(save_dir)
    path.mkdir(parents=True, exist_ok=True)
    results_filename = os.path.join(save_dir, 'results_all.csv')
    results_all_df.to_csv(results_filename, index=False)


def _update_fold_results(X, y, y_test, y_pred, test, fold, features, features_train, features_test, model):
    results = pd.DataFrame({'gene_name': y.iloc[test].index.get_level_values('gene_name').values,
                            'sample': y.iloc[test].index.get_level_values('sample').values,
                            'tissue': y.iloc[test].index.get_level_values('tissue').values,
                            'y_pred': y_pred,
                            'y_test': np.array([i for l in y_test.values.tolist() for i in l]),
                            'fold': fold})

    X_test_all = X[features].iloc[test].fillna(0)
    for feature in features:
        results[feature] = X_test_all.iloc[:, features.index(feature)].values

    if features_train != features_test:
        X_test_on_train_features = X[features_train].iloc[test].fillna(0)
        y_pred_on_train_features = model.predict_proba(X_test_on_train_features)[:, 1]
        results['y_pred_on_train_features'] = y_pred_on_train_features  
    
    return results

def _update_results_with_missing(X_missing, y_missing, y_pred_missing_iterm, features, features_train, features_test, models_all, results_all_df):
    y_pred_missing_iterm_on_train_features = y_pred_missing_iterm
    for i,m in enumerate(models_all):
        y_pred_missing_iterm[:,i] = m.predict_proba(X_missing[features_test])[:, 1]
        if features_train != features_test:
            y_pred_missing_iterm_on_train_features[:,i] = m.predict_proba(X_missing[features_train])[:, 1]

    y_pred_missing = np.mean(y_pred_missing_iterm, axis=1)
    if features_train != features_test:
        y_pred_missing_on_train_features = np.mean(y_pred_missing_iterm_on_train_features, axis=1)

    results_missing = pd.DataFrame({
                            'gene_name': y_missing.index.get_level_values('gene_name').values,
                            'sample': y_missing.index.get_level_values('sample').values,
                            'tissue': y_missing.index.get_level_values('tissue').values,
                            'y_pred': y_pred_missing,
                            'y_test': np.array([i for l in y_missing.values.tolist() for i in l]),
                            'fold': '',
    })
    for feature in features:
        results_missing[feature] = X_missing.values[:, features.index(feature)]

    if features_train != features_test:
        results_missing['y_pred_on_train_features'] = y_pred_missing_on_train_features

    results_all_df = pd.concat([results_all_df,
                                results_missing]).reset_index()
    
    return results_all_df