import sys
import pandas as pd
from rpy2.robjects import r, pandas2ri, FloatVector, ListVector, numpy2ri
import reservoirpy
import json
import matplotlib.pyplot as plt
import os
import numpy as np
from utils import generate_experiment_paths
from plot_hyperOpt import plot_hyperopt_report_v2
from hypersearch import research_v2
from hyperopt import hp

iteration = 0
def redirect_output_to_file(log_file):
    sys.stdout = open(log_file, 'a')

def reset_output():
    sys.stdout = sys.__stdout__


def convert_trials_to_dict(trials):
    """Convertit un objet Trials en dictionnaire JSON-sérialisable."""
    return [{key: (value.tolist() if isinstance(value, np.ndarray) else value) for key, value in trial.items()}
            for trial in trials.results] if hasattr(trials, "results") else []

def contains_nan_or_inf(value):
    """Vérifie si une valeur contient NaN ou Inf."""
    if isinstance(value, np.ndarray):
        return np.isnan(value).any() or np.isinf(value).any()
    elif isinstance(value, (float, np.float32, np.float64)):
        return np.isnan(value) or np.isinf(value)
    return False

def nmf_lbfgsb_hyperOpt(dataset, W_prime=None, p_prime=None, lambda_=10, gamma=100):
    """Appelle la fonction R pour effectuer l'optimisation."""
    numpy2ri.activate()
    pandas2ri.activate()
    r('library(DICEPRO)')
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
    """Objective function for HyperOpt."""
    lambda_ = kwargs.get("lambda_")
    p_prime = kwargs.get("p_prime")
    W_prime = kwargs.get("W_prime", 0)

    # Handle gamma depending on which search space is used
    if "gamma_factor" in kwargs:
        gamma_factor = kwargs["gamma_factor"]
        gamma = lambda_ * gamma_factor
    else:
        gamma_factor = None  # Not used in restrictionEspace
        gamma = kwargs.get("gamma")

    result = nmf_lbfgsb_hyperOpt(dataset, W_prime, p_prime, lambda_, gamma)
    result_dict = {name: np.array(result[i]) for i, name in enumerate(result.names)}

    if contains_nan_or_inf(result_dict.get('loss', np.array([]))) or contains_nan_or_inf(result_dict.get('constraint', np.array([]))):
        return {
            'loss': float('inf'),
            'constraint': float('inf'),
            'status': 'fail',
            'current_params': {
                'gamma_factor': gamma_factor,
                'gamma': gamma,
                'lambda_': lambda_,
                'p_prime': p_prime,
                'frobNorm': float('inf'),
                'constNorm': float('inf'),
                'c1': float('inf'),
                'c2': float('inf'),
                'objectiveValue': float('inf'),
                'penalty': float('inf')
            }
        }

    return {
        'loss': float(result_dict['objectiveValue'][0]),
        'constraint': float(result_dict['constraint'][0]),
        'status': 'OK',
        'current_params': {
            'gamma_factor': gamma_factor,
            'gamma': gamma,
            'lambda_': lambda_,
            'p_prime': p_prime,
            'frobNorm': float(result_dict['frobNorm'][0]),
            'constNorm': float(result_dict['constNorm'][0]),
            'c1': float(result_dict['c1'][0]),
            'c2': float(result_dict['c2'][0]),
            'objectiveValue': float(result_dict['objectiveValue'][0]),
            'penalty': float(result_dict['penalty'][0])
        }
    }


def custom_space():
    return {
        "lambda_": hp.loguniform("lambda_", np.log(1), np.log(1e5)),
        "gamma_factor": hp.loguniform("gamma_factor", np.log(2), np.log(1e2)),
        "p_prime": hp.loguniform("p_prime", np.log(1e-1), np.log(1)),
    }

def run_experiment(dataset, bulkName="", refName="", hp_max_evals=100, algo_select="random", output_base_dir=".", hspaceTechniqueChoose="restrictionEspace"):
    """
    Runs a hyperparameter optimization experiment and saves results in a structured directory.
    """
    pandas2ri.activate()

    # Generate necessary paths for saving experiment artifacts
    output_dir, optim_dir, data_dir, config_path, best_json_path, trials_json_path, report_dir, report_path = generate_experiment_paths(
        output_base_dir, bulkName, refName, hspaceTechniqueChoose, algo_select
    )

    hp_space = None  # Default hyperparameter space (may be customized below)

    # Choose hyperparameter search space based on the selected technique
    if hspaceTechniqueChoose == "all":
        hyperopt_config = {
            "exp": optim_dir, 
            "hp_max_evals": hp_max_evals,
            "hp_method": algo_select,
            "seed": 4,
            "hp_space": {
                "lambda_": ["loguniform", 1, 1e+8],
                "gamma": ["loguniform", 1, 1e+8],
                "p_prime": ["loguniform", 1e-6, 1],
            }
        }

    elif hspaceTechniqueChoose == "restrictionEspace": 
        hyperopt_config = {
            "exp": optim_dir,
            "hp_max_evals": hp_max_evals,
            "hp_method": algo_select,
            "seed": 4
        }
        hp_space = custom_space()  # Custom-defined hyperparameter space
    
    else:
        # Raise an error if the selected technique is invalid
        raise ValueError(f"Technique '{hspaceTechniqueChoose}' not available. Choose 'all' or 'restrictionEspace'.")

    # Save the experiment configuration as a JSON file
    with open(config_path, "w") as f:
        json.dump(hyperopt_config, f)
        
    # Run the optimization process
    best, trials = research_v2(objective, dataset, config_path, hp_space, report_dir)

    # Save all trials in JSON format
    with open(trials_json_path, "w") as f:
        json.dump(convert_trials_to_dict(trials), f, indent=4)
    
    # Plot the optimization report depending on the technique used
    if hspaceTechniqueChoose == "all":
        fig_constraint = plot_hyperopt_report_v2(optim_dir, ("gamma", "lambda_", "p_prime"), metric="constraint")
    else:  # for 'restrictionEspace'
        fig_constraint = plot_hyperopt_report_v2(optim_dir, ("gamma_factor", "lambda_", "p_prime"), metric="constraint")

    # Save the plot and display it
    fig_constraint.savefig(report_path)
    plt.show()



