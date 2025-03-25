import pandas as pd
from rpy2.robjects import r, pandas2ri, FloatVector, ListVector
import json
import matplotlib.pyplot as plt
import os
import numpy as np
import lightgbm
import seaborn as sns
from hyperopt import hp, fmin, atpe, Trials
import sys

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

def plot_results(trials, bulkName, refName, data_dir):
    """ Génère et sauvegarde les plots d'analyse de l'optimisation. """
    results = pd.DataFrame([
        {
            'loss': t['result'].get('loss', None),
            'constraint': t['result'].get('constraint', None),
            'lambda_': t['misc']['vals'].get('lambda_', [None])[0],
            'gamma': t['misc']['vals'].get('gamma', [None])[0],
            'p_prime': t['misc']['vals'].get('p_prime', [None])[0]
        }
        for t in trials.trials if t['result'].get('status') == 'ok'
    ])
    results.to_csv(f"{data_dir}/hyperopt_results_{bulkName}_{refName}.csv", index=False)
    plt.figure(figsize=(8, 5))
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(results['loss'], marker='o', linestyle='-', color='b', label='Loss')
    ax1.set_xlabel("Itérations")
    ax1.set_ylabel("Loss", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax2 = ax1.twinx()
    ax2.plot(results['constraint'], marker='s', linestyle='--', color='r', label='Contrainte')
    ax2.set_ylabel("Contrainte", color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    plt.title(f"Évolution de la Loss et de la Contrainte ({bulkName} vs {refName})")
    fig.tight_layout()
    plt.grid(True)
    plt.savefig(f"{data_dir}/loss_constraint_evolution_{bulkName}_{refName}.png")
    plt.close()
    print(f"Plots sauvegardés dans {data_dir}")

def kalidou_lbfgsb(dataset, W_prime=0, p_prime=0, lambda_=10, gamma=100):
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

def objective(dataset, bulkName, refName, config=None, **kwargs):
    """Objective function for HyperOpt with hyperparameter display."""
    global iteration
    iteration += 1
    try:
        log_file = f"dataTestNMFOptHyper/{bulkName}_{refName}/log_{bulkName}_{refName}.txt"
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        redirect_output_to_file(log_file)
        
        lambda_ = kwargs.get("lambda_")
        gamma = kwargs.get("gamma")
        p_prime = kwargs.get("p_prime")
        W_prime = kwargs.get("W_prime", 0)
        
        print(f"[Iteration {iteration}] Testing hyperparameters:")
        print(f"   - lambda_ = {lambda_:.6f}")
        print(f"   - gamma = {gamma:.6f}")
        print(f"   - p_prime = {p_prime:.6f}")
        print(f"   - W_prime = {W_prime:.6f}")
        
        result = kalidou_lbfgsb(dataset, W_prime, p_prime, lambda_, gamma)
        result_dict = {name: np.array(result[i]) for i, name in enumerate(result.names)}
        
        if contains_nan_or_inf(result_dict.get('loss', np.array([]))) or contains_nan_or_inf(result_dict.get('constraint', np.array([]))):
            print(f"❌ [Iteration {iteration}] Error: NaN or Inf detected in the results.")
            reset_output()
            return {'loss': float('inf'), 'constraint': float('inf'), 'status': 'fail'}
        
        loss = float(result_dict['loss'][0])
        constraint = float(result_dict['constraint'][0])
        
        if loss < 0:
            print(f"❌ [Iteration {iteration}] Loss is less than 0.")
            reset_output()
            return {'loss': float('inf'), 'constraint': constraint, 'status': 'fail'}
        
        if abs(constraint) > 1:
            print(f"❌ [Iteration {iteration}] Constraint has an absolute value greater than 1.")
            reset_output()
            return {'loss': loss, 'constraint': constraint, 'status': 'fail'}
        
        # 📌 Valid result
        print(f"✅ [Iteration {iteration}] Result:")
        print(f"   - loss = {loss:.6f}")
        print(f"   - constraint = {constraint:.6f}")
        print("----------------------------------------------------")
        
        reset_output()
        return {
            'loss': loss,
            'constraint': constraint,
            'status': 'ok'
        }
    
    except Exception as e:
        print(f" [Iteration {iteration}] Error in the objective function: {e}")
        reset_output()
        return {'loss': float('inf'), 'constraint': float('inf'), 'status': 'fail'}


def run_experiment(B_input, W_cb, P_cb, bulkName="", refName="", hp_max_evals=100):
    """Exécute l'expérience d'optimisation des hyperparamètres."""
    pandas2ri.activate()
    data_dir = f"dataTestNMFOptHyper/{bulkName}_{refName}"
    exp_name = f"kalidou-optim_{bulkName}_{refName}"
    
    hp_space = {
        "lambda_": hp.loguniform("lambda_", np.log(1), np.log(1e5)),  
        "gamma": hp.loguniform("gamma", np.log(1), np.log(1e5)),
        "p_prime": hp.loguniform("p_prime", np.log(0.1), np.log(1)),
    }
    if isinstance(B_input, np.ndarray):
        B_input = pd.DataFrame(B_input)
    
    if isinstance(W_cb, np.ndarray):
        W_cb = pd.DataFrame(W_cb)
    intersectionGene = B_input.index.intersection(W_cb.index)
    dataset = {'B': B_input.loc[intersectionGene], 'W': W_cb.loc[intersectionGene], 'P': P_cb}
    
    trials = Trials()
    best_results = fmin(
        fn=lambda params: objective(dataset, bulkName, refName, **params),
        space=hp_space,
        algo=atpe.suggest,
        max_evals=hp_max_evals,
        trials=trials
    )
    
    os.makedirs(data_dir, exist_ok=True)
    
    best_json_path = f"{data_dir}/best_hyperparameters_{bulkName}_{refName}.json"
    with open(best_json_path, "w") as f:
        json.dump(best_results, f, indent=4)
    
    trials_json_path = f"{data_dir}/trials_{bulkName}_{refName}.json"
    with open(trials_json_path, "w") as f:
        json.dump(convert_trials_to_dict(trials), f, indent=4)
    
    print(f"Résultats sauvegardés dans {data_dir}")
    plot_results(trials, bulkName, refName, data_dir)
