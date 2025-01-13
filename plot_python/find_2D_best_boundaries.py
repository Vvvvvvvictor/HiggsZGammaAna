import numpy as np
import pandas as pd
import uproot
import pickle as pkl
import json
import os
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from xgboost import XGBClassifier
from scipy.optimize import minimize
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import QuantileTransformer
import optuna

from pdb import set_trace
from tqdm import tqdm

# General settings
TUNING_HYPERPARAMER = True
USE_STORED_PARAMS = True
USE_STORED_MODEL = False
CATEGORIZATION = True
APPLY_BDT = True
FLOATB = True
REWEIGHT = True
TRANSFORM = True
METRIC_METHOD = 'accuracy'

# Function to calculate significance considering weights
def calculate_significance(signal, background, sig_err=None, bkg_err=None):
    ntot = signal + background
    significance = 2 * ((ntot * np.log(ntot / background)) - signal)
    if sig_err is None or bkg_err is None:
        return significance, -1
    error = (np.log(ntot/background)*sig_err)**2 + ((np.log(ntot/background) - (signal/background))*bkg_err)**2
    return significance, error

def accuracy_score(y_true, y_pred, weights, H_mass):
    H_mass_mask = (H_mass > 120) & (H_mass < 130)
    slice_thresh = np.percentile(y_pred[H_mass_mask], np.linspace(0, 100, 11))
    significance = 0
    for i in range(1, len(slice_thresh) - 1):
        mask = (y_pred >= slice_thresh[i - 1]) & (y_pred < slice_thresh[i])
        signal = np.sum(weights[(y_true == 1) & mask & H_mass_mask])
        background = np.sum(weights[(y_true == 0) & mask & H_mass_mask])
        if background > 2:
            significance += calculate_significance(signal, background)[0]
    return np.sqrt(significance)

# Example data input from user
filepath = "/eos/user/j/jiehan/root/outputs/two_jet"
bkg_names = ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
sig_names = ["ggH_M125", "VBF_M125", "ZH_M125", "WplusH_M125", "WminusH_M125", "ZH_M125", "ttH_M125"]
data_names = ["Data"]

print("Loading data...")

# Load background data
bkg_data = pd.DataFrame()
for bkg in bkg_names:
    bkg_data = pd.concat([bkg_data, uproot.open(f"{filepath}/{bkg}.root")["two_jet"].arrays(["bdt_score", "vbf_score", "weight", "H_mass", "event"], library="pd").query("H_mass > 100 & H_mass < 180")], ignore_index=True)
bkg_data["label"] = 0

# Load signal data
sig_data = pd.DataFrame()
for sig in sig_names:
    sig_data = pd.concat([sig_data, uproot.open(f"{filepath}/{sig}.root")["two_jet"].arrays(["bdt_score", "vbf_score", "weight", "H_mass","event"], library="pd").query("H_mass > 100 & H_mass < 180")], ignore_index=True)
sig_data["label"] = 1

# Load data data
data_data = pd.DataFrame()
for data in data_names:
    data_data = pd.concat([data_data, uproot.open(f"{filepath}/{data}.root")["two_jet"].arrays(["bdt_score", "vbf_score", "weight", "H_mass", "event"], library="pd").query("(H_mass > 100 & H_mass < 120) | (H_mass > 130 & H_mass < 180)")], ignore_index=True)

print("Preparing data for BDT training...")

X_train, y_train, weights_train, H_mass_train = [], [], [], []
X_test, y_test, weights_test, weights_raw_test, H_mass_test, events_test, bdt_score_test, vbf_score_test = [], [], [], [], [], [], [], []
X_val, y_val, weights_val, weights_raw_val, H_mass_val = [], [], [], [], []
X_test_data, y_test_data, weights_raw_test_data = [], [], []
bkg_H_mass_mask = (bkg_data["H_mass"] > 120) & (bkg_data["H_mass"] < 130)
sig_H_mass_mask = (sig_data["H_mass"] > 120) & (sig_data["H_mass"] < 130)
for i in range(4):
    def get_indices(data, i, mask):
        train_index = ((data["event"] % 314159 % 4 != i) & (data["event"] % 314159 % 4 != (i + 1) % 4))
        test_index = (data["event"] % 314159 % 4 == i) 
        val_index = (data["event"] % 314159 % 4 == (i + 1) % 4)
        return train_index, test_index, val_index

    bkg_train_index, bkg_test_index, bkg_val_index = get_indices(bkg_data, i, bkg_H_mass_mask)
    sig_train_index, sig_test_index, sig_val_index = get_indices(sig_data, i, sig_H_mass_mask)

    def prepare_data(bkg_index, sig_index, bkg_data, sig_data, features):
        X = np.vstack((np.concatenate((bkg_data["bdt_score"][bkg_index], sig_data["bdt_score"][sig_index])),
                       np.concatenate((bkg_data["vbf_score"][bkg_index], sig_data["vbf_score"][sig_index])))).T
        y = np.concatenate((np.zeros(len(bkg_data["bdt_score"][bkg_index])), np.ones(len(sig_data["bdt_score"][sig_index]))))
        weights = np.concatenate((bkg_data["weight"][bkg_index] / np.sum(bkg_data["weight"][bkg_index]) * sig_data["weight"][sig_index].shape[0],
                                  sig_data["weight"][sig_index] / np.mean(sig_data["weight"][sig_index])))
        return X, y, abs(weights)

    X_train_i, y_train_i, weights_train_i = prepare_data(bkg_train_index, sig_train_index, bkg_data, sig_data, ["bdt_score", "vbf_score"])
    X_test_i, y_test_i, weights_test_i = prepare_data(bkg_test_index, sig_test_index, bkg_data, sig_data, ["bdt_score", "vbf_score"])
    X_val_i, y_val_i, weights_val_i = prepare_data(bkg_val_index, sig_val_index, bkg_data, sig_data, ["bdt_score", "vbf_score"])

    X_train.append(X_train_i)
    y_train.append(y_train_i)
    weights_train.append(weights_train_i)
    H_mass_train.append(np.concatenate((bkg_data["H_mass"][bkg_train_index], sig_data["H_mass"][sig_train_index])))

    X_test.append(X_test_i)
    y_test.append(y_test_i)
    weights_test.append(weights_test_i)
    weights_raw_test.append(np.concatenate((bkg_data["weight"][bkg_test_index], sig_data["weight"][sig_test_index])))
    H_mass_test.append(np.concatenate((bkg_data["H_mass"][bkg_test_index], sig_data["H_mass"][sig_test_index])))
    events_test.append(np.concatenate((bkg_data["event"][bkg_test_index], sig_data["event"][sig_test_index])))
    bdt_score_test.append(np.concatenate((bkg_data["bdt_score"][bkg_test_index], sig_data["bdt_score"][sig_test_index])))
    vbf_score_test.append(np.concatenate((bkg_data["vbf_score"][bkg_test_index], sig_data["vbf_score"][sig_test_index])))

    X_val.append(X_val_i)
    y_val.append(y_val_i)
    weights_val.append(weights_val_i)
    weights_raw_val.append(np.concatenate((bkg_data["weight"][bkg_val_index], sig_data["weight"][sig_val_index])))
    H_mass_val.append(np.concatenate((bkg_data["H_mass"][bkg_val_index], sig_data["H_mass"][sig_val_index])))
    
    data_data_test_index = (data_data["event"] % 314159 % 4 == i)
    X_test_data.append(np.vstack((data_data["bdt_score"][data_data_test_index], data_data["vbf_score"][data_data_test_index])).T)
    y_test_data.append(np.zeros(len(data_data["bdt_score"][data_data_test_index])))
    weights_raw_test_data.append(data_data["weight"][data_data_test_index])

accuracies = []
test_set_predictions = []
test_data = []

if not os.path.exists("models"):
    os.makedirs("models")
    
def get_auc(y_true, y_pred, weights):
    fpr, tpr, _ = roc_curve(y_true, y_pred, sample_weight=weights)
    tpr, fpr = np.array(list(zip(*sorted(zip(tpr, fpr)))))
    return 1 - auc(tpr, fpr)

def plot_auc(train_y_true, train_y_pred, train_weights, val_y_true, val_y_pred, val_weights, test_y_true, test_y_pred, test_weights, fold):
    fpr_train, tpr_train, _ = roc_curve(train_y_true, train_y_pred, sample_weight=train_weights)
    fpr_val, tpr_val, _ = roc_curve(val_y_true, val_y_pred, sample_weight=val_weights)
    fpr_test, tpr_test, _ = roc_curve(test_y_true, test_y_pred, sample_weight=test_weights)
    tpr_train, fpr_train = np.array(list(zip(*sorted(zip(tpr_train, fpr_train)))))
    tpr_val, fpr_val = np.array(list(zip(*sorted(zip(tpr_val, fpr_val)))))
    tpr_test, fpr_test = np.array(list(zip(*sorted(zip(tpr_test, fpr_test)))))
    train_auc = 1 - auc(tpr_train, fpr_train)
    val_auc = 1 - auc(tpr_val, fpr_val)
    test_auc = 1 - auc(tpr_test, fpr_test)
    
    fig, ax = plt.subplots()
    ax.plot(fpr_train, tpr_train, label=f'Train (AUC = {train_auc:.3f})')
    ax.plot(fpr_val, tpr_val, label=f'Validation (AUC = {val_auc:.3f})')
    ax.plot(fpr_test, tpr_test, label=f'Test (AUC = {test_auc:.3f})')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.legend()
    plt.savefig('figs/roc_curve_fold_{}.png'.format(fold))
    
def plot_feature_importance(xgb, feature_names, fold):
    fig, ax = plt.subplots()
    ax.barh(feature_names, xgb.feature_importances_)
    ax.set_xlabel('Feature importance')
    ax.set_ylabel('Feature')
    plt.savefig('figs/feature_importance_fold_{}.png'.format(fold))

if TUNING_HYPERPARAMER:
    print("Tuning hyperparameters...")

    # Hyperparameter optimization using Optuna
    def objective(trial):
        n_estimators = trial.suggest_int('n_estimators', 50, 200)
        learning_rate = trial.suggest_float('learning_rate', 0.01, 0.3)
        max_depth = trial.suggest_int('max_depth', 3, 10)
        min_child_weight = trial.suggest_int('min_child_weight', 1, 10)
        subsample = trial.suggest_float('subsample', 0.5, 1.0)
        colsample_bytree = trial.suggest_float('colsample_bytree', 0.5, 1.0)
        
        significances_train, significances_val = [], []
        eval_aucs, train_aucs = [], []
        
        for i in range(4):         
            xgb = XGBClassifier(
                n_estimators=n_estimators,
                learning_rate=learning_rate,
                max_depth=max_depth,
                min_child_weight=min_child_weight,
                subsample=subsample,
                colsample_bytree=colsample_bytree,
                objective='binary:logistic',
                eval_metric=['auc', 'logloss'],
                early_stopping_rounds=20
            )
                               
            xgb.fit(X_train[i], y_train[i], sample_weight=weights_train[i], eval_set=[(X_val[i], y_val[i])], verbose=False)
            y_pred = xgb.predict_proba(X_val[i])[:, 1]
            significance_val = accuracy_score(y_val[i], y_pred, weights_raw_val[i], H_mass_val[i])
            significances_val.append(significance_val)
            significance_train = accuracy_score(y_train[i], xgb.predict_proba(X_train[i])[:, 1], weights_train[i], H_mass_train[i])
            significances_train.append(significance_train)
            
            eval_aucs.append(get_auc(y_val[i], y_pred, weights_raw_val[i]))
            train_aucs.append(get_auc(y_train[i], xgb.predict_proba(X_train[i])[:, 1], weights_train[i]))
            
        optuna_metrics = {
            'accuracy': np.sqrt(np.mean(train_aucs)*(2*np.mean(eval_aucs) - np.mean(train_aucs))),
            'significance': np.sum(np.array(significances_val)**2)**0.5
        }
        return optuna_metrics[METRIC_METHOD]

    study = optuna.create_study(direction='maximize')
    study.optimize(objective, n_trials=40)

    best_params = study.best_params
    print(f"Best parameters: {best_params}")

    with open("models/best_params.json", "w") as f:
        json.dump(best_params, f)

if os.path.exists("models/best_params.json") and USE_STORED_PARAMS:
    best_params = json.load(open("models/best_params.json", "r"))

print("Training models...")

xgbs, tfers = [], []
 # Train and evaluate the model with the best hyperparameters
for i in range(4):
    if os.path.exists(f"models/xgb_models_{i}.pkl") and USE_STORED_MODEL:
        xgb = pkl.load(open(f"models/xgb_models_{i}.pkl", "rb"))
    else:
        # Initialize and train model with best hyperparameters
        if USE_STORED_PARAMS:
            xgb = XGBClassifier(**best_params, objective='binary:logistic', eval_metric=['auc', 'logloss'], early_stopping_rounds=20)
        else:
            xgb = XGBClassifier(objective='binary:logistic', eval_metric=['auc', 'logloss'], early_stopping_rounds=20)
        xgb.fit(X_train[i], y_train[i], sample_weight=weights_train[i], eval_set=[(X_val[i], y_val[i])], verbose=False)
        pkl.dump(xgb, open(f"models/xgb_models_{i}.pkl", "wb"))
    
    xgbs.append(xgb)
    plot_feature_importance(xgb, ['Dijet', 'VBF'], i)
    plot_auc(y_train[i], xgb.predict_proba(X_train[i])[:, 1], weights_train[i], y_val[i], xgb.predict_proba(X_val[i])[:, 1], weights_val[i], y_test[i], xgb.predict_proba(X_test[i])[:, 1], weights_raw_test[i], i)
    
    # Transform the score
    tfer = QuantileTransformer(output_distribution='uniform', random_state=0, n_quantiles=1000)
    tfer.fit(xgb.predict_proba(X_test[i])[y_test[i] == 1].reshape(-1, 1))
    tfers.append(tfer)
    
    # Predict and evaluate
    y_pred = xgb.predict_proba(X_test[i])[:, 1]
    accuracy = accuracy_score(y_test[i], y_pred, weights_raw_test[i], H_mass_test[i])
    accuracies.append(accuracy)
    print(f"Fold accuracy: {accuracy}")
    
    # Store the test set predictions and corresponding weights
    test_set_predictions.append(pd.DataFrame({
        'bdt_score': X_test[i][:, 0],
        'vbf_score': X_test[i][:, 1],
        'label': y_test[i],
        'H_mass': H_mass_test[i],
        'weight': weights_raw_test[i],
        'score': y_pred,
        'score_t': tfer.transform(y_pred.reshape(-1, 1)).flatten()
    }))
    
    score_data = xgb.predict_proba(X_test_data[i])[:, 1]
    test_data.append(pd.DataFrame({
        'bdt_score': X_test_data[i][:, 0],
        'vbf_score': X_test_data[i][:, 1],
        'weight': weights_raw_test_data[i],
        'score': score_data,
        'score_t': tfer.transform(score_data.reshape(-1, 1)).flatten()
    }))

print(f"Total significance: {np.sum(np.array(accuracies)**2)**0.5}")

# Combine all test set predictions
test_set_predictions_df = pd.concat(test_set_predictions)
test_data = pd.concat(test_data)


def calculate_weights(data, bl, br, data_err=None):
    weight = data[bl:br].sum()
    if data_err is None:
        return weight, -1
    weight_err = np.sqrt(data_err[bl:br].sum())
    return weight, weight_err

# Function to optimize thresholds for maximum significance
def fit(sig, sig_err, bkg, bkg_err, bl, br, nbin, minN=2, early_stop=-1, floatB=False, pbar=False, num=1, den=1):
    if nbin == 1:
        if floatB: return [], 0, 0

        signal, signal_err = calculate_weights(sig, bl, br, sig_err)
        background, background_err = calculate_weights(bkg, bl, br, bkg_err)
        
        if REWEIGHT:
            rw_num, _ = calculate_weights(num, bl, br)
            rw_den, _ = calculate_weights(den, bl, br)
            if rw_num < 10 or rw_den < 10: return -1, -1, -1
            background = background * rw_num / rw_den
        
        if background < minN:
            return -1, -1, -1
        significance, error = calculate_significance(signal, background, signal_err, background_err)
        return [bl], significance, error

    elif nbin > 1:
        L = int(np.ceil(np.log2(nbin)))
        N2 = 2**(L-1)
        N1 = nbin - N2
        
        b_opt, s_opt, e_opt, stop = -1, -1, -1, 0
        
        for b in (tqdm(np.arange(bl, br)) if pbar else np.arange(bl, br)):
            b1, s1, e1 = fit(sig, sig_err, bkg, bkg_err, bl, b, N1, minN, floatB=floatB, num=num, den=den)
            if b1 == -1:
                continue
            b2, s2, e2 = fit(sig, sig_err, bkg, bkg_err, b, br, N2, minN, num=num, den=den)
            if b2 == -1:
                break
            if s1 + s2 > s_opt:
                s_opt = s1 + s2
                b_opt = sorted(list(set(b1 + [b] + b2)))
                e_opt = (e1*s1) + (e2*s2)
            else:
                stop += 1
                if stop == early_stop:
                    break

        return b_opt, s_opt, e_opt

# Find the best n_cat by optimizing the thresholds for maximum significance
n_scan = 100
if CATEGORIZATION:
    best_n_cat = 4
    best_significance, best_significance_error = 0, 0
    for n_cat in range(2, 11):  # Try different numbers of categories from 2 to 10
        sig_query_str = 'H_mass > 120 & H_mass < 130 & label == 1'
        bkg_query_str = 'H_mass > 120 & H_mass < 130 & label == 0'
        if REWEIGHT:
            den_query_str = '(H_mass > 100 & H_mass < 120) | (H_mass > 130 & H_mass < 180)'
            optimized_thresholds, total_significance, significance_error = fit(
                np.histogram(test_set_predictions_df.query(sig_query_str)['score_t'], bins=n_scan, range=(0, 1), weights=test_set_predictions_df.query(sig_query_str)['weight'])[0],
                np.histogram(test_set_predictions_df.query(sig_query_str)['score_t'], bins=n_scan, range=(0, 1), weights=test_set_predictions_df.query(sig_query_str)['weight']**2)[0],
                np.histogram(test_set_predictions_df.query(bkg_query_str)['score_t'], bins=n_scan, range=(0, 1), weights=test_set_predictions_df.query(bkg_query_str)['weight'])[0],
                np.histogram(test_set_predictions_df.query(bkg_query_str)['score_t'], bins=n_scan, range=(0, 1), weights=test_set_predictions_df.query(bkg_query_str)['weight']**2)[0],
                0, n_scan, (n_cat + 1 if FLOATB else n_cat),
                pbar=True, floatB=FLOATB,
                num=np.histogram(test_data['score_t'], bins=n_scan, range=(0, 1), weights=test_data['weight'])[0],
                den=np.histogram(test_set_predictions_df.query(den_query_str)['score_t'], bins=n_scan, range=(0, 1), weights=test_set_predictions_df.query(den_query_str)['weight'])[0]
            )
        else:
            optimized_thresholds, total_significance, significance_error = fit(
                np.histogram(test_set_predictions_df.query(sig_query_str)['score_t'], bins=n_scan, range=(0, 1))[0],
                np.histogram(test_set_predictions_df.query(bkg_query_str)['score_t'], bins=n_scan, range=(0, 1))[0],
                0, n_scan, (n_cat + 1 if FLOATB else n_cat),
                pbar=True, floatB=FLOATB
            )
        optimized_thresholds.append(n_scan)
        optimized_thresholds = list(map(round, np.array(optimized_thresholds)/n_scan, [3]*len(optimized_thresholds)))
        total_significance = np.sqrt(total_significance)
        significance_error = np.sqrt(significance_error/total_significance)
        
        print(f"Number of categories: {n_cat}, total significance: {total_significance} +/- {significance_error}, previous best: {best_significance}, thresholds: {optimized_thresholds}")
        if (total_significance - best_significance) / np.sqrt(significance_error**2 + best_significance_error**2) > 0.4750: 
            # Require at least p((s-s_best)>0) > 68.26%
        # if total_significance > best_significance * 1.01:  # Require at least 1% improvement
            best_significance = total_significance
            best_significance_error = significance_error
            best_n_cat = n_cat
            best_thresholds = optimized_thresholds
        if n_cat > best_n_cat + 1:  # Stop if the significance does not improve for 1 consecutive categories
            break

    print(f"Best number of categories: {best_n_cat}")
    print(f"Best total significance: {best_significance}")
    with open("models/best_thresholds.txt", "w") as f:
        f.write(" ".join(map(str, best_thresholds)))
else:
    with open("models/best_thresholds.txt", "r") as f:
        best_thresholds = list(map(float, f.readline().split()))
        best_n_cat = len(best_thresholds) - 1

# Assign category labels to data based on the best thresholds
cat_labels = np.digitize(test_set_predictions_df[f"score{'_t' if TRANSFORM else ''}"], best_thresholds) - 1
test_set_predictions_df['category'] = cat_labels

print("Significance in each category:")
for i in range(best_n_cat):
    selection = f"label == 1 & H_mass > 120 & H_mass < 130 & score{'_t' if TRANSFORM else ''} >= {best_thresholds[i]} & score{'_t' if TRANSFORM else ''} < {best_thresholds[i + 1]}"
    signal = np.sum(test_set_predictions_df.query(selection)['weight'])
    background = np.sum(test_set_predictions_df.query(selection.replace("label == 1", "label == 0"))['weight'])
    signal_err = np.sum(test_set_predictions_df.query(selection)['weight']**2)**0.5
    background_err = np.sum(test_set_predictions_df.query(selection.replace("label == 1", "label == 0"))['weight']**2)**0.5
    if REWEIGHT:
        num = test_data.query(selection.replace("label == 1 & H_mass > 120 & H_mass < 130 & ", ""))
        num = np.histogram(num['score_t'], bins=n_scan, range=(0, 1), weights=num['weight'])[0]
        den = test_set_predictions_df.query(f"label == 0 & ((H_mass > 100 & H_mass < 120) | (H_mass > 130 & H_mass < 180)) & score{'_t' if TRANSFORM else ''} >= {best_thresholds[i]} & score{'_t' if TRANSFORM else ''} < {best_thresholds[i + 1]}")
        den = np.histogram(den['score_t'], bins=n_scan, range=(0, 1), weights=den['weight'])[0]
        rw_num, _ = calculate_weights(num, int(best_thresholds[i]*n_scan), int(best_thresholds[i + 1]*n_scan))
        rw_den, _ = calculate_weights(den, int(best_thresholds[i]*n_scan), int(best_thresholds[i + 1]*n_scan))
        background = background * rw_num / rw_den
        background_err = background_err * (rw_num / rw_den)
    significance, error = calculate_significance(signal, background, signal_err, background_err)
    print(f"Category {i}: {np.sqrt(significance)} +/- {np.sqrt(error/significance)}, signal: {signal} +/- {signal_err}, background: {background} +/- {background_err}")

with uproot.recreate("output_with_categories.root") as f:
    f["test"] = {
        'bdt_score': test_set_predictions_df['bdt_score'].astype(np.float32),
        'vbf_score': test_set_predictions_df['vbf_score'].astype(np.float32),
        'label': test_set_predictions_df['label'].astype(np.int32),
        'H_mass': test_set_predictions_df['H_mass'].astype(np.float32),
        'weight': test_set_predictions_df['weight'].astype(np.float32),
        'score': test_set_predictions_df['score'].astype(np.float32),
        'score_t': test_set_predictions_df['score_t'].astype(np.float32),
        'category': test_set_predictions_df['category'].astype(np.int32)
    }
    
if APPLY_BDT:
    def apply_bdt(xgbs, name, boundaries, input_dir, output_dir):
        data = uproot.open(f"{input_dir}/{name}.root")["two_jet"].arrays(library="pd").query("H_mass > 100 & H_mass < 180")
        data.rename(columns={"bdt_score_t": "two_jet_score_t", "bdt_score": "two_jet_score"}, inplace=True)
        data_o = pd.DataFrame()
        for i in range(4):
            data_cat= data.query(f"event % 314159 % 4 == {i}")
            X = np.vstack((data_cat["two_jet_score"], data_cat["vbf_score"])).T
            data_cat["bdt_score"] = xgbs[i].predict_proba(X)[:, 1]
            data_cat["bdt_score_t"] = tfers[i].transform(data_cat["bdt_score"].values.reshape(-1, 1)).flatten()
            data_o = pd.concat([data_o, data_cat], ignore_index=True)
        data_o["category"] = np.digitize(data_o["bdt_score_t"], boundaries) - 1
        with uproot.recreate(f"{output_dir}/{name}.root") as f:
            f["two_jet"] = data_o
        
    input_dir = "/eos/user/j/jiehan/root/outputs/two_jet"
    output_dir = "/eos/user/j/jiehan/root/outputs/two_jet_2D_categories"
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
    for bkg in bkg_names:
        print(f"Applying BDT to {bkg}...")
        apply_bdt(xgbs, bkg, best_thresholds, input_dir, output_dir)
    for sig in sig_names:
        print(f"Applying BDT to {sig}...")
        apply_bdt(xgbs, sig, best_thresholds, input_dir, output_dir)
    for data in data_names:
        print(f"Applying BDT to {data}...")
        apply_bdt(xgbs, data, best_thresholds, input_dir, output_dir)
    
# plot signal and background score distributions
fig = plt.figure()
selection = 'H_mass > 120 & H_mass < 130'
var = f"score{'_t' if TRANSFORM else ''}"
signal_hist, bins = np.histogram(test_set_predictions_df.query(f'label == 1 & {selection}')[var], bins=100, range=(0, 1), weights=test_set_predictions_df.query(f'label == 1 & {selection}')['weight'])
signal_hist_err = np.histogram(test_set_predictions_df.query(f'label == 1 & {selection}')[var], bins=bins, range=(0, 1), weights=test_set_predictions_df.query(f'label == 1 & {selection}')['weight']**2)[0]**0.5 / np.sum(signal_hist)
print("signal yield: ", np.sum(signal_hist))
signal_hist = signal_hist / np.sum(signal_hist)
background_hist, bins = np.histogram(test_set_predictions_df.query(f'label == 0 & {selection}')[var], bins=bins, range=(0, 1), weights=test_set_predictions_df.query(f'label == 0 & {selection}')['weight'])
background_hist_err = np.histogram(test_set_predictions_df.query(f'label == 0 & {selection}')[var], bins=bins, range=(0, 1), weights=test_set_predictions_df.query(f'label == 0 & {selection}')['weight']**2)[0]**0.5 / np.sum(background_hist)
print("background yield: ", np.sum(background_hist))
background_hist = background_hist / np.sum(background_hist)

pos = (bins[:-1] + bins[1:]) / 2
plt.errorbar(pos, signal_hist, yerr=signal_hist_err, fmt='o', label='Signal', color='b')
plt.errorbar(pos, background_hist, yerr=background_hist_err, fmt='o', label='Background', color='r')

for i in range(len(best_thresholds)):
    plt.axvline(x=best_thresholds[i], color='k', linestyle='--')

plt.xlabel('Categorization score')
plt.ylabel('Events')
plt.legend()
plt.savefig('figs/score_distribution.png')

fig = plt.figure()
selection = 'H_mass < 120 | H_mass > 130'
var = f"score{'_t' if TRANSFORM else ''}"
background_hist, bins = np.histogram(test_set_predictions_df.query(f'label == 0 & {selection}')[var], bins=100, range=(0, 1), weights=test_set_predictions_df.query(f'label == 0 & {selection}')['weight'])
background_hist_err = np.histogram(test_set_predictions_df.query(f'label == 0 & {selection}')[var], bins=bins, range=(0, 1), weights=test_set_predictions_df.query(f'label == 0 & {selection}')['weight']**2)[0]**0.5 / np.sum(background_hist)
print("background yield: ", np.sum(background_hist))
background_hist = background_hist / np.sum(background_hist)
data_hist, bins = np.histogram(test_data[var], bins=bins, range=(0, 1), weights=test_data['weight'])
data_hist_err = np.histogram(test_data[var], bins=bins, range=(0, 1), weights=test_data['weight']**2)[0]**0.5 / np.sum(data_hist) 
print("data yield: ", np.sum(data_hist))
data_hist = data_hist / np.sum(data_hist)

pos = (bins[:-1] + bins[1:]) / 2
plt.errorbar(pos, background_hist, yerr=background_hist_err, fmt='o', label='Background', color='r')
plt.errorbar(pos, data_hist, yerr=data_hist_err, fmt='o', label='Data', color='g')

for i in range(len(best_thresholds)):
    plt.axvline(x=best_thresholds[i], color='k', linestyle='--')

plt.xlabel('Categorization score')
plt.ylabel('Events')
plt.legend()
plt.savefig('figs/score_distribution_data.png')