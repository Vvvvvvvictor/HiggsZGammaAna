import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted

class WeightedQuantileTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, n_quantiles=1000, output_distribution='uniform', random_state=None, subsample=int(1e9)):
        self.n_quantiles = n_quantiles
        if output_distribution != 'uniform':
            raise ValueError("WeightedQuantileTransformer currently only supports 'uniform' output_distribution.")
        self.output_distribution = output_distribution
        self.random_state = random_state # For potential subsampling, not implemented yet
        self.subsample = subsample # Not implemented yet
        # self.quantiles_ stores the target quantiles, e.g., [0, 0.1, ..., 1.0] for uniform output
        self.quantiles_ = np.linspace(0, 1, self.n_quantiles)


    def fit(self, X, y=None, sample_weight=None):
        X = check_array(X, ensure_2d=False, dtype="numeric", copy=True)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        if X.shape[1] > 1:
            raise ValueError("WeightedQuantileTransformer only supports 1D input.")

        scores = X[:, 0]
        
        if sample_weight is None:
            sample_weight = np.ones_like(scores)
        else:
            sample_weight = np.asarray(sample_weight).flatten() # Ensure 1D

        if np.sum(sample_weight) <= 0: # All weights are zero or negative, or no samples
            min_score, max_score = (np.min(scores) if len(scores) > 0 else 0), (np.max(scores) if len(scores) > 0 else 0)
            if min_score == max_score: 
                 min_score -= 1e-7 # Avoid singular interval
                 max_score += 1e-7
            self.references_ = np.linspace(min_score, max_score, self.n_quantiles)
            if len(np.unique(self.references_)) < 2 and self.n_quantiles > 1 : 
                 self.references_[0] = min_score -1e-7
                 self.references_[-1] = max_score + 1e-7 if self.n_quantiles > 1 else max_score 
                 if self.n_quantiles == 1: self.references_[0] = min_score 

            self.is_fitted_ = True
            return self

        finite_mask = np.isfinite(scores) & np.isfinite(sample_weight) & (sample_weight > 0)
        
        if not np.any(finite_mask): # No valid data points
            min_s, max_s = (np.min(X[:,0]) if X.shape[0] > 0 else 0), (np.max(X[:,0]) if X.shape[0] > 0 else 0)
            if min_s == max_s: 
                min_s -=1e-7
                max_s += 1e-7
            self.references_ = np.linspace(min_s, max_s, self.n_quantiles)
            if len(np.unique(self.references_)) < 2 and self.n_quantiles > 1 :
                 self.references_[0] = min_s -1e-7
                 self.references_[-1] = max_s + 1e-7 if self.n_quantiles > 1 else max_s
                 if self.n_quantiles == 1: self.references_[0] = min_s
            self.is_fitted_ = True
            return self

        scores_finite = scores[finite_mask]
        sample_weight_finite = sample_weight[finite_mask]

        sorted_indices = np.argsort(scores_finite)
        _scores = scores_finite[sorted_indices]
        _weights = sample_weight_finite[sorted_indices]
        
        cumulative_weights = np.cumsum(_weights)
        total_weight = cumulative_weights[-1]

        if total_weight <= 0: 
            min_s_val, max_s_val = np.min(_scores) if len(_scores) > 0 else 0, np.max(_scores) if len(_scores) > 0 else 0
            if min_s_val == max_s_val and len(_scores) > 0:
                min_s_val -= 1e-7
                max_s_val += 1e-7
            elif len(_scores) == 0:
                min_s_val, max_s_val = 0, 1 
            self.references_ = np.linspace(min_s_val, max_s_val, self.n_quantiles)
            self.is_fitted_ = True
            return self

        ecdf_y = cumulative_weights / total_weight
        self.references_ = np.interp(self.quantiles_, ecdf_y, _scores)
        
        if len(np.unique(self.references_)) < 2 and self.n_quantiles > 1:
            if len(np.unique(_scores)) > 1 : 
                 pass 
            else: 
                 ref_val = self.references_[0]
                 self.references_[0] = ref_val - 1e-7
                 self.references_[-1] = ref_val + 1e-7
                 if self.n_quantiles > 2: 
                     self.references_ = np.linspace(ref_val - 1e-7, ref_val + 1e-7, self.n_quantiles)

        self.is_fitted_ = True
        return self

    def transform(self, X):
        check_is_fitted(self, "is_fitted_")
        X = check_array(X, ensure_2d=False, dtype="numeric", copy=True)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        if X.shape[1] > 1:
            raise ValueError("WeightedQuantileTransformer only supports 1D input.")

        input_scores = X[:, 0]
        transformed_scores = np.interp(input_scores, self.references_, self.quantiles_, left=0., right=1.)
        
        return transformed_scores.reshape(-1, 1)
