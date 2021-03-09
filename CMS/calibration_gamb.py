#! /usr/bin/env python3

from argparse import ArgumentParser
import csv
from datetime import datetime
import logging
import json
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

import optuna
import cmsmodel
import hatmodel_ug_gamb
import hatmodel_ss
import hatmodel_drc_fusc
import hatmodel_drc_quan

builders = {
   "hatmodel_ug_gamb": hatmodel_ug_gamb.build_model,
   "hatmodel_ss": hatmodel_ss.build_model,
   "hatmodel_drc_fusc": hatmodel_drc_fusc.build_model,
   "hatmodel_drc_quan": hatmodel_drc_quan.build_model
}

import clr
# Normally I would do this down in if __name__ == "__main__", but we need the path to compartments[.exe] for the clr.AddReference call below.
parser = ArgumentParser()
# The default value here will work if the .NET assembly "compartments" is in the PYTHONPATH.
# If you are using the pycms docker container, this will be the case. Note that the default value
# doesn't have ".exe" at the end of it.
parser.add_argument("-c", "--compartments", default="bin/compartments", help="Specify full path to compartments.exe")
parser.add_argument("-b", "--database", type=str, default="sqlite:///hat-studies.db", help="Study database name ['hat-studies.db']")
parser.add_argument("-d", "--data", type=Path, default=Path("HAT_ug_Tbg_calibr.csv"), help="Target case data for calibration")
parser.add_argument("-m", "--model", type=str, default="hatmodel_ug_gamb", help=f"model name - [{', '.join(list(builders.keys()))}]")
parser.add_argument("-n", "--name", type=str, default="default", help="Study name ['default']")
parser.add_argument("-p", "--png", action="store_true", help="Save output to a .png file")
parser.add_argument("-r", "--repetitions", type=int, default=3, help="Number of CMS solver repetitions per trial")
parser.add_argument("-s", "--surveillance", type=float, default=0.65, help="Estimated surveillance rate")
parser.add_argument("-t", "--trials", type=int, default=10, help="Number of trials in study [10]")
parser.add_argument("-v", "--verbosity", type=str, default="INFO", help="{DEBUG|INFO|WARN|ERROR|CRITICAL}")
parser.add_argument("-w", "--working", type=Path, default=None, help="Working directory path")
parser.add_argument("-y", "--years", type=int, default=45, help="Number of years in simulation [45]")

args = parser.parse_args()
args.working = args.working if args.working else Path.cwd().absolute() / "studies" / f"{datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M')}-{args.name}"
logger = logging.getLogger(__name__)

logger.info(f"Loading CMS from {args.compartments}")
clr.AddReference(args.compartments)
from compartments.emodl import EmodlLoader
from compartments import Configuration as cfg
from compartments.emod.utils import SolverFactory as solvers


def configure_logging(level: str = "INFO", directory: Path = None):

    levels = {
        "DEBUG": logging.DEBUG,         # 10
        "INFO": logging.INFO,           # 20
        "WARN": logging.WARN,           # 30
        "WARNING": logging.WARNING,     # 30
        "ERROR": logging.ERROR,         # 40
        "CRITICAL": logging.CRITICAL,   # 50
        "FATAL": logging.FATAL          # 50
    }

    level = levels[level] if isinstance(level, str) else int(level)

    # logger = logging.getLogger()    # did this at module scope above
    logger.setLevel(level)
    console = logging.StreamHandler()
    console.setLevel(level)
    formatter = logging.Formatter("%(levelname)s:%(message)s")
    console.setFormatter(formatter)
    logger.addHandler(console)
    logfile = (directory / f"{__file__}-{datetime.now():%Y%m%d-%H%M%S}.log").absolute()
    print(f"Logging to {logfile}")
    disk = logging.FileHandler(logfile)
    disk.setLevel(level)
    disk.setFormatter(formatter)
    logger.addHandler(disk)

    return


args.working.mkdir(exist_ok=True)
configure_logging(args.verbosity, args.working)

df = pd.read_csv(args.data)
case_data = np.array(df.groupby("Year")["New_HAT_cases"].sum(), dtype=np.uint32)


def objective(trial):

    beta_h = trial.suggest_uniform('beta_h', 0.001, 0.1)

    config = {
        "solver": "TAU",
        # "prng_seed": 20201025,  # use this for stochasticity: datetime.now().microsecond
        "prng_seed": datetime.now().microsecond,
        "tau-leaping": {
            "epsilon": 0.001,
            "Nc": 10,
            "Multiple": 10.0,
            "SSARuns": 100
        }
    }
    cfg.CurrentConfiguration = cfg.ConfigurationFromString(json.dumps(config))

    model_name = args.model.lower()
    if model_name in builders:
        model_description = builders[model_name](beta_h, **{
            "human-susceptible":2_030-546-189,
            "human-infectious-one":546,
            "human-infectious-two":189,
            "tsetse-susceptible":2_030*6.56*(1-0.006818574),
            "tsetse-infectious":2_030*6.56*0.006818574,
            "reservoir-susceptible":2_030*0.24*(1-0.2847),
            "reservoir-infectious":2_030*0.24*0.2847,
            "non-reservoir-hosts":2_030/10})
    else:
        raise RuntimeError(f"Model{model_name}is not a known model.")

    model_info = load_model(model_description)

    # Create an SSA solver - could just specify "SSA" or "TAU" here. Run for 3650 time units, record 3650 intermediate states of the system.
    args.repetitions = 3
    args.years = 45
    sim_duration = args.years*365 + 1   # 5 years + 1 more day
    num_samples = sim_duration
    t_start = datetime.now()
    solver = solvers.CreateSolver(config["solver"], model_info, args.repetitions, sim_duration, num_samples)
    t_create = datetime.now()
    logger.info(f"{t_create - t_start} for creating the CMS solver.")
    solver.Solve()  # Run the solver
    t_solve = datetime.now()
    logger.info(f"{t_solve - t_create} for solving the model {args.repetitions} times ({sim_duration} simulated days).")

    datafile = args.working / f"trajectories-{trial.number:03}.csv"
    save_trajectories_to_file(solver, datafile)

    save_plots(solver, trial, args.working)

    # extract relevant trajectories
    trajectories = []
    data = solver.GetTrajectoryData()   # Retrieve the recorded data (all observables, see build_model())
    max_sample = num_samples - 1
    extra_data = max_sample - (args.years * 365)
    print(f"Ignoring {extra_data} additional samples") if extra_data else None
    year_indices = np.array(range(args.years+1))*365
    for index, label in enumerate(solver.GetTrajectoryLabels()):
        if label.startswith("human-infection-cumulative"):
            trajectory = np.array(list(data[index]), dtype=np.uint32)
            trajectory = trajectory[year_indices]

            # transform cumulative counts to incidence (t - t-1)
            trajectory[1:] -= trajectory[0:-1]

            trajectories.append(trajectory)

    # aggregate counts by year
    trajectories = np.array(trajectories)
    # take mean of all trajectories
    mean = trajectories[0, :] if trajectories.shape[0] == 1 else np.mean(trajectories, axis=0)
    # fitness is sum(squares) last N years
    num_years_data = len(case_data)
    fitness = case_data - (mean[-num_years_data:] * args.surveillance)
    fitness *= fitness
    score = np.sum(fitness)

    return score


def build_model(beta_h = 0.0076, **kwargs):

    # See hatmodel.md

    model = cmsmodel.CmsModel("hat")

    # The "odd" construction of `kwargs[...] if ... in kwargs else 0` allows you to selectively specify
    # some initial populations in the call to build_model() and those values will be used to initialize
    # the population of those species. If you do not specify a value, the initial population will be 0.
    species = [
        {"name": "human-susceptible",          "population": kwargs["human-susceptible"]          if "human-susceptible"          in kwargs else 0, "observe":True},
        {"name": "human-exposed",              "population": kwargs["human-exposed"]              if "human-exposed"              in kwargs else 0, "observe":True},
        {"name": "human-infectious-one",       "population": kwargs["human-infectious-one"]       if "human-infectious-one"       in kwargs else 0, "observe":True},
        {"name": "human-infectious-two",       "population": kwargs["human-infectious-two"]       if "human-infectious-two"       in kwargs else 0, "observe":True},
        {"name": "human-recovered",            "population": kwargs["human-recovered"]            if "human-recovered"            in kwargs else 0, "observe":True},
        {"name": "human-infection-cumulative", "population": kwargs["human-infection-cumulative"] if "human-infection-cumulative" in kwargs else 0, "observe":True},
        {"name": "human-dead",                 "population": kwargs["human-dead"]                 if "human-dead"                 in kwargs else 0, "observe":True},

        {"name": "tsetse-susceptible",     "population": kwargs["tsetse-susceptible"]     if "tsetse-susceptible"     in kwargs else 0, "observe":True},
        {"name": "tsetse-exposed",         "population": kwargs["tsetse-exposed"]         if "tsetse-exposed"         in kwargs else 0, "observe":True},
        {"name": "tsetse-infectious",      "population": kwargs["tsetse-infectious"]      if "tsetse-infectious"      in kwargs else 0, "observe":True},
        {"name": "tsetse-non-susceptible", "population": kwargs["tsetse-non-susceptible"] if "tsetse-non-susceptible" in kwargs else 0, "observe":True},

        {"name": "reservoir-susceptible", "population": kwargs["reservoir-susceptible"] if "reservoir-susceptible" in kwargs else 0, "observe":True},
        {"name": "reservoir-exposed",     "population": kwargs["reservoir-exposed"]     if "reservoir-exposed"     in kwargs else 0, "observe":True},
        {"name": "reservoir-infectious",  "population": kwargs["reservoir-infectious"]  if "reservoir-infectious"  in kwargs else 0, "observe":True},
        {"name": "reservoir-recovered",   "population": kwargs["reservoir-recovered"]   if "reservoir-recovered"   in kwargs else 0, "observe":True}
    ]

    def _add_species(name: str, population: int, observe: bool):
        model.add_species(name, population, observe)

    for specie in species:
        _add_species(**specie)

    parameters = [
        {"name":"sigma-h",          "value":0.083333333},    # incubation rate (human E->I1)
        {"name":"phi-h",            "value":0.001901141},    # progression from human I1->I2
        {"name":"omega-h",          "value":0.005479452},     # rate of progression from I2-> death
        {"name":"beta-v",           "value":0.212},     # beta for tsetse fly infection from infectious human or reservoir animal
        {"name":"p-human-feed",     "value":0.089},      # probability of human feed. Changed to 0.2 in South Sudan, just assumed
        #{"name":"p-reservoir-feed", "value":0.153},      # probability of reservoir host feed
        {"name":"sigma-v",          "value":0.06},     # incubation rate (tsetse E->I)
        {"name":"mu-v",             "value":0.03846154},      # tsetse fly mortality rate
        {"name":"treatment-one",    "value":0.1608275},      # probability a human is treated in stage 1
        {"name":"treatment-two",    "value":0.5384225},      # probability a human is treated in stage 2
        # _not_ from the paper referenced above
        {"name":"p-feed",           "value":1.0/3},     # probability of feeding in a given 24 hours
        {"name":"beta-h",           "value":beta_h},       # beta for human infection by infectious tsetse fly. Changed to 0.26 for South Sudan, just assumed
        {"name":"beta-r",           "value":0.1345},       # beta for reservoir infection by infectious tsetse fly
        {"name":"phi-r",            "value":0.14285714},    # reservoir incubation rate, swine
        {"name":"omega-r-nt",       "value":1/182.5},     # reservoir recovery rate with treatment (assume no recovery otherwise); assume treatment q6 months. FIT
        {"name":"omega-r-t",        "value":1/91.25},     # reservoir recovery rate without treatment
        {"name":"treatment-reservoir",          "value":0},     # probability reservoir is treated (assume lifelong infection otherwise)
        {"name":"mu-h",             "value":0.000053},  # human mortality rate
        {"name":"mu-r",             "value":0.001369863},    # reservoir mortality rate
        {"name":"wane_immunity",             "value":1/50}    # same for animals and humans
    ]

    def _add_parameter(name: str, value: float):
        model.add_parameter(name, value)

    for parameter in parameters:
        _add_parameter(**parameter)

    # Convenience functions:
    model.add_function("human-population",           "(+ human-susceptible human-exposed human-infectious-one human-infectious-two human-recovered)")
    model.add_function("reservoir-population",       "(+ reservoir-susceptible reservoir-exposed reservoir-infectious reservoir-recovered)")
    model.add_function("tsetse-population",          "(+ tsetse-susceptible tsetse-exposed tsetse-infectious tsetse-non-susceptible)")
    model.add_function("tsetse-human-ratio",         "(/ tsetse-population human-population)")
    model.add_function("tsetse-reservoir-ratio",     "(/ tsetse-population reservoir-population)")

    # Reactions/transitions

    # human S->E->I1->I2->R
    model.add_reaction("human-infection",              ["human-susceptible"],    ["human-exposed", "human-infection-cumulative"], "(* human-susceptible beta-h p-feed p-human-feed tsetse-human-ratio (/ tsetse-infectious tsetse-population))")
    model.add_reaction("human-exposed-infectious",     ["human-exposed"],        ["human-infectious-one"],                        "(* sigma-h human-exposed)")
    model.add_reaction("human-infectious-progression", ["human-infectious-one"], ["human-infectious-two"],                        "(* phi-h human-infectious-one (- 1 treatment-one))")
    model.add_reaction("human-recovery-passive",       ["human-infectious-two"], ["human-recovered"],                             "(* omega-h human-infectious-two treatment-two)")
    model.add_reaction("human-recovery-active",        ["human-infectious-one"], ["human-recovered"],                             "(* phi-h human-infectious-one treatment-one)")
    model.add_reaction("human-death",                  ["human-infectious-two"], ["human-dead"],                                  "(* omega-h human-infectious-two (- 1 treatment-two))")

    # Tsetse S->E->I->R, NS
    model.add_function("infectious-feed", "(+ (* p-human-feed (/ human-infectious-one human-population)) (* p-reservoir-feed (/ reservoir-infectious reservoir-population)))")  # probability a given feed is infectious and transmission occurs
    # model.add_function("infectious-feed", "(* beta-v p-feed (+ (* p-human-feed (/ human-infectious-one human-population)) (* p-reservoir-feed (/ reservoir-infectious reservoir-population))))")#probablility a given feed is infectious and transmission occurs
    model.add_reaction("feed-and-infected",      ["tsetse-susceptible"], ["tsetse-exposed"], "(* p-feed beta-v infectious-feed tsetse-susceptible)")    # probability a feed happens in first 24 hours, is infectious, and transmission occurs
    # model.add_reaction("become-non-susceptible", ["tsetse-susceptible"], ["tsetse-non-susceptible"], "(+ (- 1 p-feed)(* p-feed (- 1 infectious-feed)))")#probability tsetse doesn't feed in first 24 hours, or does but transmission doesn't occur
    model.add_reaction("become-non-susceptible", ["tsetse-susceptible"], ["tsetse-non-susceptible"], "(* (+ (- 1 p-feed)(* p-feed (- 1 infectious-feed))) tsetse-susceptible)")     # probability tsetse doesn't feed in first 24 hours, or does but transmission doesn't occur
    model.add_reaction("tsetse-progress-to-infectious", ["tsetse-exposed"], ["tsetse-infectious"], "(* sigma-v tsetse-exposed)")

    # reservoir S->E->I->R
    model.add_reaction("reservoir-infection",          ["reservoir-susceptible"], ["reservoir-exposed"],    "(* reservoir-susceptible beta-r p-feed p-reservoir-feed tsetse-reservoir-ratio (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious", ["reservoir-exposed"],     ["reservoir-infectious"], "(* phi-r reservoir-exposed)")
    model.add_reaction("reservoir-recovery",           ["reservoir-infectious"],  ["reservoir-recovered"],  "(* omega-r treatment-reservoir reservoir-infectious)")

    # vital dynamics: stable population - recycle deaths directly into births (susceptible)

    # model.add_reaction("human-susceptible-death-birth",    ["human-susceptible"],   ["human-susceptible"], "(* mu-h human-susceptible)")
    model.add_reaction("human-exposed-death-birth",        ["human-exposed"],        ["human-susceptible"], "(* mu-h human-exposed)")
    model.add_reaction("human-infectious-one-death-birth", ["human-infectious-one"], ["human-susceptible"], "(* mu-h human-infectious-one)")
    model.add_reaction("human-infectious-two-death-birth", ["human-infectious-two"], ["human-susceptible"], "(* mu-h human-infectious-two)")
    model.add_reaction("human-recovered-death-birth",      ["human-recovered"],      ["human-susceptible"], "(* mu-h human-recovered)")

    # model.add_reaction("vector-susceptible-death-birth",     ["tsetse-susceptible"],     ["tsetse-susceptible"], "(* mu-v tsetse-susceptible)")
    model.add_reaction("vector-exposed-death-birth",         ["tsetse-exposed"],         ["tsetse-susceptible"], "(* mu-v tsetse-exposed)")
    model.add_reaction("vector-infectious-death-birth",      ["tsetse-infectious"],      ["tsetse-susceptible"], "(* mu-v tsetse-infectious)")
    model.add_reaction("vector-non-susceptible-death-birth", ["tsetse-non-susceptible"], ["tsetse-susceptible"], "(* mu-v tsetse-non-susceptible)")

    # model.add_reaction("reservoir-susceptible-death-birth", ["reservoir-susceptible"], ["reservoir-susceptible"], "(* mu-r reservoir-susceptible)")
    model.add_reaction("reservoir-exposed-death-birth",     ["reservoir-exposed"],     ["reservoir-susceptible"], "(* mu-r reservoir-exposed)")
    model.add_reaction("reservoir-infectious-death-birth",  ["reservoir-infectious"],  ["reservoir-susceptible"], "(* mu-r reservoir-infectious)")
    model.add_reaction("reservoir-recovered-death-birth",   ["reservoir-recovered"],   ["reservoir-susceptible"], "(* mu-r reservoir-recovered)")

    return model


def load_model(model):

    model_info = EmodlLoader.LoadEMODLModel(str(model))

    return model_info


def save_trajectories_to_file(solver, datafile):
    logger.info(f"Saving trajectories to {datafile}")
    data = solver.GetTrajectoryData()
    with datafile.open("w") as file:
        writer = csv.writer(file)
        for index, label in enumerate(solver.GetTrajectoryLabels()):
            row = [label]
            row.extend(data[index])     # add actual data to row
            writer.writerow(row)        # write row to file
    return


def save_plots(solver, trial, working_directory):
    logger.info("Saving plots")
    data = solver.GetTrajectoryData()
    for fig, population in enumerate(["human", "tsetse", "reservoir"]):
        plt.figure(fig+trial.number*3)  # starts a new plot (3 = number of plots per trial)
        for index, label in enumerate(solver.GetTrajectoryLabels()):
            # only include this row/label if it is for the current population
            if label.startswith(population):
                # Each "row" of data is one trajectory of one observable, the label is observable_name{run#}
                color = ["red", "darkorange", "gold", "green", "blue", "indigo", "violet"][index % 7]
                style = ["solid", "dashed", "dashdot", "dotted"][(index // 7) % 4]
                marker = [None, "o", "^", "1", "+", "x"][(index // 28) % 6]
                plt.plot([float(value) for value in data[index]], label=str(label), color=color, linestyle=style, marker=marker)
        plt.legend()
        filename = working_directory / f"{population}-{trial.number:03}.png"
        logger.info(f"Saving plot to '{filename}'")
        plt.savefig(filename)

    return


def log_parameters(arguments):

    logger.info(f"CMS binary:            {arguments.compartments}")
    logger.info(f"model:.................{arguments.model}")
    logger.info(f"working directory:     {arguments.working}")
    logger.info(f"database: .............{arguments.database}")
    logger.info(f"study name:            {arguments.name}")
    logger.info(f"number of trials: .....{arguments.trials}")
    logger.info(f"repetitions per trial: {arguments.repetitions}")
    logger.info(f"simulation years: .....{arguments.years}")
    logger.info(f"calibration data:      {arguments.data}")
    logger.info(f"surveillance rate .....{arguments.surveillance}")

    return


if __name__ == "__main__":

    log_parameters(args)

    logger.info(f"optuna.create_study(storage='{args.database}', study_name='{args.name}', direction='minimize', load_if_exists=True)")
    study = optuna.create_study(storage=args.database, study_name=args.name, direction="minimize", load_if_exists=True)
    logger.info(f"study.optimize(objective, n_trials={args.trials})")
    study.optimize(objective, n_trials=args.trials)
    logger.info(f"best params = {study.best_params}")
    logger.info(f"best trial = {study.best_trial}")
    logger.info(f"best_value = {study.best_value}")
