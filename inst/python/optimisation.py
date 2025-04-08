import sys
import pandas as pd
from rpy2.robjects import r, pandas2ri, FloatVector, ListVector, numpy2ri
import reservoirpy
import json
import matplotlib.pyplot as plt
import os
import numpy as np
from utils import generate_experiment_paths

iteration = 0
def redirect_output_to_file(log_file):
    sys.stdout = open(log_file, 'a')

def reset_output():
    sys.stdout = sys.__stdout__


def convert_trials_to_dict(trials):
    """Convertit un objet Trials en dictionnaire JSON-sérialisable."""
    return [{key: (value.tolist() if isinstance(value, np.ndarray) else value)
             for key, value in trial.items()}
            for trial in trials.results] if hasattr(trials, "results") else []

def contains_nan_or_inf(value):
    """Vérifie si une valeur contient NaN ou Inf."""
    if isinstance(value, np.ndarray):
        return np.isnan(value).any() or np.isinf(value).any()
    elif isinstance(value, (float, np.float32, np.float64)):
        return np.isnan(value) or np.isinf(value)
    return False

def kalidou_lbfgsb(dataset, W_prime=None, p_prime=None, lambda_=10, gamma=100):
    """Appelle la fonction R pour effectuer l'optimisation."""
    numpy2ri.activate()
    pandas2ri.activate()
    r('library(DICEPRO)')
    print(type(dataset['B']))
    B_r = pandas2ri.py2rpy(pd.DataFrame(dataset['B'])) if isinstance(dataset['B'], np.ndarray) else dataset['B']
    W_cb_r = pandas2ri.py2rpy(pd.DataFrame(dataset['W'])) if isinstance(dataset['W'], np.ndarray) else dataset['W']
    P_cb_r = pandas2ri.py2rpy(pd.DataFrame(dataset['P'])) if isinstance(dataset['P'], np.ndarray) else dataset['P']
    r_dataset = ListVector({
        'B': B_r,
        'W_cb': W_cb_r,
        'P_cb': P_cb_r
    })
    W_prime_r = numpy2ri.py2rpy(W_prime) if W_prime is not None else None
    p_prime_r = numpy2ri.py2rpy(p_prime) if p_prime is not None else None
    r_lambda = FloatVector([lambda_])
    r_gamma = FloatVector([gamma])
    r_func = r['nmf_lbfgsb']
    result = r_func(r_dataset, W_prime, p_prime, r_lambda, r_gamma)
    return result

def objective(dataset, config=None, **kwargs):
    """Fonction objectif pour HyperOpt."""
    lambda_ = kwargs.get("lambda_")
    gamma = kwargs.get("gamma")
    p_prime = kwargs.get("p_prime")
    W_prime = kwargs.get("W_prime", 0)
    
    result = kalidou_lbfgsb(dataset, W_prime, p_prime, lambda_, gamma)
    result_dict = {name: np.array(result[i]) for i, name in enumerate(result.names)}
    if contains_nan_or_inf(result_dict.get('loss', np.array([]))) or contains_nan_or_inf(result_dict.get('constraint', np.array([]))):
        return {'loss': float('inf'), 'constraint': float('inf'), 'status': 'fail'}
    return {
        'loss': float(result_dict['loss'][0]),
        'constraint': float(result_dict['constraint'][0]),
        'status': 'OK'
    }

def run_experiment(dataset, bulkName="", refName="", hp_max_evals=100, algo_select="tpe", output_base_dir="."):
    """Exécute l'expérience d'optimisation des hyperparamètres et sauvegarde les résultats dans un dossier organisé."""
    print(dataset)
    pandas2ri.activate()
    output_dir, optim_dir, data_dir, config_path, best_json_path, trials_json_path, report_path = generate_experiment_paths(output_base_dir, bulkName, refName)
    
    hyperopt_config = {
        "exp": optim_dir,
        "hp_max_evals": hp_max_evals,
        "hp_method": algo_select,
        "seed": 4,
        "hp_space": {
            "lambda_": ["loguniform", 1, 1e+8],
            "gamma": ["loguniform", 1, 1e+5],
            "p_prime": ["loguniform", 1e-1, 1],
        },
    }
    with open(config_path, "w") as f:
        json.dump(hyperopt_config, f)
    best_results, trials = reservoirpy.hyper.research(objective, dataset, config_path, ".")
    with open(best_json_path, "w") as f:
        json.dump(best_results, f, indent=4)
    with open(trials_json_path, "w") as f:
        json.dump(convert_trials_to_dict(trials), f, indent=4)
    fig_constraint = reservoirpy.hyper.plot_hyperopt_report(optim_dir, ("lambda_", "gamma", "p_prime"), metric="constraint")
    fig_constraint.savefig(report_path)
    plt.show()




