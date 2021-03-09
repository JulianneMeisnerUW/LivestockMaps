#! /usr/bin/env python3

from argparse import ArgumentParser
import csv
from datetime import datetime
import json
from matplotlib import pyplot as plt
from pathlib import Path

import clr
import cmsmodel


def main():

    # Most of this is ignored because we pass explicit values into solvers.CreateSolver() below. Really, only "prng_seed" is active.
    config = {
        "solver": "SSA",
        "runs": 11,
        "duration": 180,
        "samples": 180,
        "prng_seed": 20201025   # use this for stocasticity: datetime.now().microsecond
    }
    cfg.CurrentConfiguration = cfg.ConfigurationFromString(json.dumps(config))

    # Selectively populate a few of the species. For more fine grained control, update the code in build_model().
    model_description = build_model(**{
        "human-susceptible":13_191-1595-553,
        "human-infectious-one":1595,
        "human-infectious-two":553,
        "tsetse-susceptible":13_191*6.56*(1-0.006818574), #2 looks pretty good
        "tsetse-infectious":13_191*6.56*0.006818574,
        "reservoir-susceptible":13_191*0.06*(1-0.2847),
        "reservoir-infectious":13_191*0.036*0.2847,
        "non-reservoir-hosts":13_191/10})
    model_info = load_model(model_description, cleanup=False)

    # Create an SSA solver - could just specify "SSA" or "TAU" here. Run for 3650 time units, record 3650 intermediate states of the system.
    solver = solvers.CreateSolver(config["solver"], model_info, 10, 365.0*40, 365*40)
    solver.Solve()  # Run the solver
    data = solver.GetTrajectoryData()   # Retrieve the recorded data (all observables, see build_model())

    # could include the csv file name in command line arguments
    with Path("trajectories_drc_fusc_int_25tt_25itc.csv").open("w") as file:
        writer = csv.writer(file)
        for index, label in enumerate(solver.GetTrajectoryLabels()):
            row = [label.replace("{0}", "")]    # compartments are suffixed with trajectory #, but there is only one trajectory, so clean it up
            row.extend(data[index])             # add actual data to row
            writer.writerow(row)                # write row to file

    # iterate over the three populations of interest
    for fig, population in enumerate(["human", "tsetse", "reservoir"]):
        plt.figure(fig)     # starts a new plot
        for index, label in enumerate(solver.GetTrajectoryLabels()):
            # only include this row/label if it is for the current population
            if label.startswith(population):
                # Each "row" of data is one trajectory of one observable, the label is observable_name{run#}
                color = ["red", "darkorange", "gold", "green", "blue", "indigo", "violet"][index % 7]
                style = ["solid", "dashed", "dashdot", "dotted"][(index // 7) % 4]
                marker = [None, "o", "^", "1", "+", "x"][(index // 28) % 6]
                plt.plot([float(value) for value in data[index]], label=str(label), color=color, linestyle=style, marker=marker)
        plt.legend()
        #if not args.png:
        #    plt.show()
        #else:
        #must indent the next three lines if I restore the ifelse
        # use population name for output
        print(f"Saving plot to '{population}.png'")
        plt.savefig(f"{population}.png")

    return


def build_model(beta_h=0.009680082507907695, **kwargs):

    # See hatmodel.md

    model = cmsmodel.CmsModel("hat")

    # The "odd" construction of `kwargs[...] if ... in kwargs else 0` allows you to selectively specify
    # some initial populations in the call to build_model() and those values will be used to initialize
    # the population of those species. If you do not specify a value, the initial population will be 0.
    species = [
        {"name":"human-susceptible",          "population": kwargs["human-susceptible"]          if "human-susceptible"          in kwargs else 0, "observe":True},
        {"name":"human-exposed",              "population": kwargs["human-exposed"]              if "human-exposed"              in kwargs else 0, "observe":True},
        {"name":"human-infectious-one",       "population": kwargs["human-infectious-one"]       if "human-infectious-one"       in kwargs else 0, "observe":True},
        {"name":"human-infectious-two",       "population": kwargs["human-infectious-two"]       if "human-infectious-two"       in kwargs else 0, "observe":True},
        {"name":"human-recovered",            "population": kwargs["human-recovered"]            if "human-recovered"            in kwargs else 0, "observe":True},
        {"name":"human-infection-cumulative", "population": kwargs["human-infection-cumulative"] if "human-infection-cumulative" in kwargs else 0, "observe":True},
        {"name":"human-dead",                 "population": kwargs["human-dead"]                 if "human-dead"                 in kwargs else 0, "observe":True},

        {"name":"tsetse-susceptible",     "population": kwargs["tsetse-susceptible"]     if "tsetse-susceptible"     in kwargs else 0, "observe":True},
        {"name":"tsetse-exposed",         "population": kwargs["tsetse-exposed"]         if "tsetse-exposed"         in kwargs else 0, "observe":True},
        {"name":"tsetse-infectious",      "population": kwargs["tsetse-infectious"]      if "tsetse-infectious"      in kwargs else 0, "observe":True},
        {"name":"tsetse-non-susceptible", "population": kwargs["tsetse-non-susceptible"] if "tsetse-non-susceptible" in kwargs else 0, "observe":True},

        {"name":"reservoir-susceptible", "population": kwargs["reservoir-susceptible"] if "reservoir-susceptible" in kwargs else 0, "observe":True},
        {"name":"reservoir-exposed",     "population": kwargs["reservoir-exposed"]     if "reservoir-exposed"     in kwargs else 0, "observe":True},
        {"name":"reservoir-infectious",  "population": kwargs["reservoir-infectious"]  if "reservoir-infectious"  in kwargs else 0, "observe":True},
        {"name":"reservoir-recovered",   "population": kwargs["reservoir-recovered"]   if "reservoir-recovered"   in kwargs else 0, "observe":True},

        {"name":"non-reservoir-hosts",   "population": kwargs["non-reservoir-hosts"]   if "non-reservoir-hosts"   in kwargs else 0, "observe":True}
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
        {"name":"p-human-feed",     "value":0.42},      # probability of human feed. Changed to 0.2 in South Sudan, just assumed
        #{"name":"p-reservoir-feed", "value":0.153},      # probability of reservoir host feed
        {"name":"sigma-v",          "value":0.056},     # incubation rate (tsetse E->I)
        #{"name":"mu-v",             "value":0.03846154},      # tsetse fly mortality rate
        {"name":"p-survives-feed",  "value":0.89},
        {"name":"treatment-one",    "value":0.25},      # probability a human is treated in stage 1
        {"name":"treatment-two",    "value":0.25},      # probability a human is treated in stage 2
        # _not_ from the paper referenced above
        {"name":"p-feed",           "value":1.0/3},     # probability of feeding in a given 24 hours
        {"name":"beta-h",           "value":beta_h},       # beta for human infection by infectious tsetse fly. Changed to 0.26 for South Sudan, just assumed
        {"name":"beta-r",           "value":0.1345},       # beta for reservoir infection by infectious tsetse fly
        {"name":"phi-r",            "value":0.14285714},    # reservoir incubation rate, swine
        {"name":"omega-r-nt",       "value":1/182.5},     # reservoir recovery rate with treatment (assume no recovery otherwise); assume treatment q6 months. FIT
        {"name":"omega-r-t",        "value":1/91.25},     # reservoir recovery rate without treatment
        {"name":"treatment-reservoir",          "value":0.25},     # probability reservoir is treated (assume lifelong infection otherwise)
        {"name":"mu-h",             "value":0.000051},  # human mortality rate
        {"name":"mu-r",             "value":0.001369863},    # reservoir mortality rate
        {"name":"wane_immunity",    "value":1/50},    # same for animals and humans
        {"name":"p-itc",            "value":0.25}
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
    model.add_function("p-reservoir-feed",           "(* (- 1 p-human-feed) (/ reservoir-population (+ reservoir-population non-reservoir-hosts)))")
    model.add_function("mu-v",                       "(/(-(ln (* (- 1 p-itc) p-reservoir-feed p-survives-feed)))3) ")

    # Reactions/transitions

    # human S->E->I1->I2->R
    model.add_reaction("human-infection",              ["human-susceptible"],    ["human-exposed", "human-infection-cumulative"], "(* human-susceptible beta-h p-feed p-human-feed tsetse-human-ratio (/ tsetse-infectious tsetse-population))")
    model.add_reaction("human-exposed-infectious",     ["human-exposed"],        ["human-infectious-one"],                        "(* sigma-h human-exposed)")
    model.add_reaction("human-infectious-progression", ["human-infectious-one"], ["human-infectious-two"],                        "(* phi-h human-infectious-one (- 1 treatment-one))")
    model.add_reaction("human-recovery-passive",       ["human-infectious-two"], ["human-recovered"],                             "(* omega-h human-infectious-two treatment-two)")
    model.add_reaction("human-recovery-active",        ["human-infectious-one"], ["human-recovered"],                             "(* phi-h human-infectious-one treatment-one)")
    model.add_reaction("human-waning-immunity",        ["human-recovered"], ["human-susceptible"],                                "(* wane_immunity human-recovered)")
    model.add_reaction("human-death",                  ["human-infectious-two"], ["human-dead"],                                  "(* omega-h human-infectious-two (- 1 treatment-two))")

    # Tsetse S->E->I->R, NS
    model.add_function("infectious-feed", "(+ (* p-human-feed (/ human-infectious-one human-population)) (* p-reservoir-feed (/ reservoir-infectious reservoir-population)))")#probablility a given feed is infectious and transmission occurs
    #model.add_function("infectious-feed", "(* beta-v p-feed (+ (* p-human-feed (/ human-infectious-one human-population)) (* p-reservoir-feed (/ reservoir-infectious reservoir-population))))")#probablility a given feed is infectious and transmission occurs
    model.add_reaction("feed-and-infected",      ["tsetse-susceptible"], ["tsetse-exposed"], "(* p-feed beta-v infectious-feed tsetse-susceptible)") #probabilty a feed happens in first 24 hours, is infectious, and transmission occurs
    #model.add_reaction("become-non-susceptible", ["tsetse-susceptible"], ["tsetse-non-susceptible"], "(+ (- 1 p-feed)(* p-feed (- 1 infectious-feed)))")#probability tsetse doesn't feed in first 24 hours, or does but transmission doesn't occur
    model.add_reaction("become-non-susceptible", ["tsetse-susceptible"], ["tsetse-non-susceptible"], "(* (+ (- 1 p-feed)(* p-feed (- 1 infectious-feed))) tsetse-susceptible)")#probability tsetse doesn't feed in first 24 hours, or does but transmission doesn't occur
    model.add_reaction("tsetse-progress-to-infectious", ["tsetse-exposed"], ["tsetse-infectious"], "(* sigma-v tsetse-exposed)")

    # reservoir S->E->I->R
    model.add_reaction("reservoir-infection",          ["reservoir-susceptible"], ["reservoir-exposed"],    "(* reservoir-susceptible beta-r p-feed p-reservoir-feed tsetse-reservoir-ratio (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious", ["reservoir-exposed"],     ["reservoir-infectious"], "(* phi-r reservoir-exposed)")
    model.add_reaction("reservoir-recovery",           ["reservoir-infectious"],  ["reservoir-recovered"],  "(+ (* omega-r-t treatment-reservoir reservoir-infectious) (* omega-r-nt (- 1 treatment-reservoir) reservoir-infectious))")
    model.add_reaction("reservoir-waning-immunity",    ["reservoir-recovered"],   ["reservoir-susceptible"], "(* wane_immunity reservoir-recovered)")

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


def load_model(model, cleanup=True):

    model_info = EmodlLoader.LoadEMODLModel(str(model))

    return model_info


if __name__ == "__main__":

    parser = ArgumentParser()
    # The default value here will work if the .NET assembly "compartments" is in the PYTHONPATH.
    # If you are using the pycms docker container, this will be the case. Note that the default value
    # doesn't have ".exe" at the end of it.
    parser.add_argument("-c", "--compartments", default="bin/compartments",
                        help="Specify full path to compartments.exe")
    parser.add_argument("-p", "--png", action="store_true", help="Save output to a .png file")

    args = parser.parse_args()

    clr.AddReference(args.compartments)
    from compartments.emodl import EmodlLoader
    from compartments import Configuration as cfg
    from compartments.emod.utils import SolverFactory as solvers

    main()
