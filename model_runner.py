# Information/Social Epidemiology
#
# SPEC collaborative
#
# 2021

from batchrunner_local import BatchRunnerMP
from multiprocessing import freeze_support
from isepi import ISEpiModel
from isepi import Stage
from isepi import MobilityType
import pandas as pd
import json
import sys

num_procs = int(sys.argv[1])
file_params = sys.argv[2]

# Read JSON file
with open(file_params) as f:
  data = json.load(f)

print(f"Scenario description")
print(f"Description: { data['description'] }")
print(f"Prepared by: { data['prepared-by'] }")
print(f"Date: { data['date'] }")
print("")
print("Attempting to configure model from file...")

model_params = {
    #"prop_initial_infected": data["model"]["epidemiology"]["prop_initial_infected"],
    #"avg_recovery_time": data["model"]["epidemiology"]["avg_recovery_time"],
    "epidemiology": data["model"]["epidemiology"]
}

var_params = {"dummy": range(25,50,25)}

num_iterations = data["ensemble"]["runs"]
num_steps = data["ensemble"]["days"]

if __name__ == "__main__":
    freeze_support()

    batch_run = BatchRunnerMP(
        ISEpiModel,
        nr_processes=num_procs,
        fixed_parameters=model_params,
        variable_parameters=var_params,
        iterations=num_iterations,
        max_steps=num_steps*ISEpiModel.dwell_15_day,
        model_reporters={
                    "Step": compute_stepno,
                    "Susceptible": compute_susceptible,
                    "Incubating": compute_incubating,
                    "Asymptomatic": compute_asymptomatic,
                    "SymptQuarantined": compute_symptdetected,
                    "AsymptQuarantined": compute_asymptdetected,
                    "Severe": compute_severe,
                    "Recovered": compute_recovered,
                    "Deceased": compute_deceased,
                    "CummulPrivValue": compute_cumul_private_value,
                    "CummulPublValue": compute_cumul_public_value,
                    "CummulTestCost": compute_cumul_testing_cost,
                    "Rt": compute_eff_reprod_number,
                    "Employed": compute_employed,
                    "Unemployed": compute_unemployed
                },
        display_progress=True)

    print("Parametrization complete:")
    print("")
    print(json.dumps(data, indent=3))
    print("")
    print(f"Executing an ensemble of size {num_iterations} using {num_steps} steps with {num_procs} machine cores...")
    cm_runs = batch_run.run_all()

    print("")
    print("Saving results to file...")

    ldfs = []
    i = 0

    for cm in cm_runs.values():
        cm["Iteration"] = i
        ldfs.append(cm)
        i = i + 1

    file_out = data["output"]["prefix"]

    dfs = pd.concat(ldfs)
    dfs.to_csv(file_out + ".csv")

    print("Simulation completed without errors.")
