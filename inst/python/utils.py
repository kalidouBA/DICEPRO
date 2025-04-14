import os
import json

def generate_experiment_paths(output_dir, bulkName, refName):
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, f"{bulkName}_{refName}")
    os.makedirs(data_dir, exist_ok=True)
    optim_dir = os.path.join(data_dir, "optim")
    os.makedirs(optim_dir, exist_ok=True)
    report_dir = os.path.join(optim_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    exp_name = "optim"
    config_path = os.path.join(optim_dir, f"{exp_name}.config.json")
    best_json_path = os.path.join(optim_dir, f"best_hyperparameters.json")
    trials_json_path = os.path.join(optim_dir, f"trials.json")
    report_path = os.path.join(report_dir, f"hyperopt_report_constraint.png")
    return output_dir, optim_dir, data_dir, config_path, best_json_path, trials_json_path, report_dir, report_path
