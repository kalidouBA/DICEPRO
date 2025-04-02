import sys
import pandas as pd
from rpy2.robjects import r, pandas2ri, FloatVector, ListVector
import reservoirpy
import json
import matplotlib.pyplot as plt
import os
import numpy as np

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
    pandas2ri.activate()
    r('library(DICEPRO)')
    r_dataset = ListVector({
        'B': pandas2ri.py2rpy(dataset['B']),
        'W_cb': pandas2ri.py2rpy(dataset['W']),
        'P_cb': pandas2ri.py2rpy(dataset['P'])
    })
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

def run_experiment(dataset, bulkName="", refName="", hp_max_evals=100, algo_select="tpe"):
    """Exécute l'expérience d'optimisation des hyperparamètres."""
    pandas2ri.activate()
    data_dir = f"dataTestNMFOptHyper/{bulkName}_{refName}"
    exp_name = f"optim_{bulkName}_{refName}"
    hyperopt_config = {
        "exp": exp_name,
        "hp_max_evals": hp_max_evals,
        "hp_method": algo_select,
        "seed": 4,
        "hp_space": {
            "lambda_": ["loguniform", 1, 1e+5],
            "gamma": ["loguniform", 1, 1e+6],
            "p_prime": ["loguniform", 1e-1, 1],
        },
    }

    os.makedirs(data_dir, exist_ok=True)
    config_path = f"{data_dir}/{exp_name}.config.json"
    with open(config_path, "w") as f:
        json.dump(hyperopt_config, f)
        
    best_results, trials = reservoirpy.hyper.research(objective, dataset, config_path, ".")
    best_json_path = f"{data_dir}/best_hyperparameters_{bulkName}_{refName}.json"
    with open(best_json_path, "w") as f:
        json.dump(best_results, f, indent=4)
    trials_json_path = f"{data_dir}/trials_{bulkName}_{refName}.json"
    with open(trials_json_path, "w") as f:
        json.dump(convert_trials_to_dict(trials), f, indent=4)
    fig_constraint = reservoirpy.hyper.plot_hyperopt_report(exp_name, ("lambda_", "gamma", "p_prime"), metric="constraint")
    fig_constraint.savefig(f"{data_dir}/hyperopt_report_constraint_{bulkName}_{refName}.png")
    plt.show()
    
