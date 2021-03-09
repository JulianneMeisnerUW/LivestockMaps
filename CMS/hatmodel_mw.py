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
    model_description = build_model(**{ #uncomment this for running calibration_mw.py
    #model_description = build_model(**{
                "human-susceptible":939-26-70,
                "human-infectious-one":26,
                "human-infectious-two": 70,
                "tsetse-susceptible":939*6.56*(1-0.022),
                #"tsetse-susceptible":8_312, #10_807; 837*13*(1-0.006818574), same tsetse-human ratio as gambiense areas, but p_h is different/from literature
                "tsetse-infectious":939*6.56*0.022, #74; 837*13*0.006818574,
                #"tsetse-infectious":57, #74; 837*13*0.006818574,
                "reservoir-sd-susceptible":939*0.06*(1-0.013547841),
                "reservoir-sd-infectious": 939*0.06*0.013547841,
                "reservoir-rd-susceptible": 939*0.16*(1-0.061),
                "reservoir-rd-infectious": 939*0.16*0.061,
                "reservoir-sw-susceptible": 939*(0.06/2)*(1-0.017688), #later do the ratio as fitted, currently just setting to 1/2 domestic #s
                "reservoir-sw-infectious": 939*(0.06/2)*0.017688, #later do the ratio as fitted, currently just setting to 1/2 domestic #s
                "reservoir-rw-susceptible": 939*(0.16/2)*(1-0.021505376),#later do the ratio as fitted, currently just setting to 1/2 domestic #s
                "reservoir-rw-infectious": 939*(0.16/2)*0.021505376,#later do the ratio as fitted, currently just setting to 1/2 domestic #s
                "non-reservoir-hosts": 939/10})
    model_info = load_model(model_description, cleanup=False)

    # Create an SSA solver - could just specify "SSA" or "TAU" here. Run for 3650 time units, record 3650 intermediate states of the system.
    solver = solvers.CreateSolver(config["solver"], model_info, 10, 365.0*40, 365*40)
    solver.Solve()  # Run the solver
    data = solver.GetTrajectoryData()   # Retrieve the recorded data (all observables, see build_model())

    # could include the csv file name in command line arguments
    with Path("trajectories_mw.csv").open("w") as file:
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


def build_model(beta_h=0.000912069, **kwargs):

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

        {"name":"reservoir-sd-susceptible", "population": kwargs["reservoir-sd-susceptible"] if "reservoir-sd-susceptible" in kwargs else 0, "observe":True},
        {"name":"reservoir-sd-exposed",     "population": kwargs["reservoir-sd-exposed"]     if "reservoir-sd-exposed"     in kwargs else 0, "observe":True},
        {"name":"reservoir-sd-infectious",  "population": kwargs["reservoir-sd-infectious"]  if "reservoir-sd-infectious"  in kwargs else 0, "observe":True},
        {"name":"reservoir-sd-recovered",   "population": kwargs["reservoir-sd-recovered"]   if "reservoir-sd-recovered"   in kwargs else 0, "observe":True},

        {"name":"reservoir-rd-susceptible", "population": kwargs["reservoir-rd-susceptible"] if "reservoir-rd-susceptible" in kwargs else 0, "observe":True},
        {"name":"reservoir-rd-exposed",     "population": kwargs["reservoir-rd-exposed"]     if "reservoir-rd-exposed"     in kwargs else 0, "observe":True},
        {"name":"reservoir-rd-infectious",  "population": kwargs["reservoir-rd-infectious"]  if "reservoir-rd-infectious"  in kwargs else 0, "observe":True},
        {"name":"reservoir-rd-recovered",   "population": kwargs["reservoir-rd-recovered"]   if "reservoir-rd-recovered"   in kwargs else 0, "observe":True},

        {"name":"reservoir-sw-susceptible", "population": kwargs["reservoir-sw-susceptible"] if "reservoir-sw-susceptible" in kwargs else 0, "observe":True},
        {"name":"reservoir-sw-exposed",     "population": kwargs["reservoir-sw-exposed"]     if "reservoir-sw-exposed"     in kwargs else 0, "observe":True},
        {"name":"reservoir-sw-infectious",  "population": kwargs["reservoir-sw-infectious"]  if "reservoir-sw-infectious"  in kwargs else 0, "observe":True},
        {"name":"reservoir-sw-recovered",   "population": kwargs["reservoir-sw-recovered"]   if "reservoir-sw-recovered"   in kwargs else 0, "observe":True},

        {"name":"reservoir-rw-susceptible", "population": kwargs["reservoir-rw-susceptible"] if "reservoir-rw-susceptible" in kwargs else 0, "observe":True},
        {"name":"reservoir-rw-exposed",     "population": kwargs["reservoir-rw-exposed"]     if "reservoir-rw-exposed"     in kwargs else 0, "observe":True},
        {"name":"reservoir-rw-infectious",  "population": kwargs["reservoir-rw-infectious"]  if "reservoir-rw-infectious"  in kwargs else 0, "observe":True},
        {"name":"reservoir-rw-recovered",   "population": kwargs["reservoir-rw-recovered"]   if "reservoir-rw-recovered"   in kwargs else 0, "observe":True},

        {"name":"non-reservoir-hosts",   "population": kwargs["non-reservoir-hosts"]   if "non-reservoir-hosts"   in kwargs else 0, "observe":True}
    ]

    def _add_species(name: str, population: int, observe: bool):
        model.add_species(name, population, observe)

    for specie in species:
        _add_species(**specie)

    parameters = [
        {"name":"sigma-h",          "value":0.083333333},    # incubation rate (human E->I1)
        {"name":"phi-h",            "value":1/30},    # progression from human I1->I2
        {"name":"omega-h",          "value":1/(82-30)},     # rate of progression from I2-> death
        {"name":"beta-v",           "value":0.258},     # beta for tsetse fly infection from infectious human or reservoir animal
        {"name":"p-human-feed",     "value":0.05},      # probability of human feed
        #{"name":"p-reservoir-feed", "value":0.781},      # probability of reservoir host feed
        {"name":"sigma-v",          "value":0.06},     # incubation rate (tsetse E->I)
        {"name":"mu-v",             "value":0.03448276},      # tsetse fly mortality rate
        {"name":"treatment-one",    "value":0},      # probability a human is treated in stage 1
        {"name":"treatment-two",    "value":0.083333333},      # probability a human is treated in stage 2
        {"name":"p-feed",           "value":1.0/3},     # probability of feeding in a given 24 hours
        {"name":"beta-h",           "value":beta_h},       # beta for human infection by infectious tsetse fly
        {"name":"beta-r",           "value":0.1345},       # beta for reservoir infection by infectious tsetse fly; FIT
        {"name":"phi-r-s",            "value":0.14285714},    # reservoir incubation rate, swine
        {"name":"phi-r-r",            "value":0.083333333},    # reservoir incubation rate, ruminants
        {"name":"omega-r-nt-s",       "value":1/182.5},     # reservoir recovery rate with treatment (assume no recovery otherwise); assume treatment q6 months. FIT
        {"name":"omega-r-nt-r",       "value":1/225},     # reservoir recovery rate with treatment (assume no recovery otherwise); assume treatment q6 months. FIT
        {"name":"omega-r-t",          "value":1/91.25},     # reservoir recovery rate without treatment (assume no recovery otherwise); assume treatment q6 months. FIT
        {"name":"treatment-reservoir-d",          "value":0.5},     # probability domestic reservoir is treated (assume lifelong infection otherwise)
        {"name":"treatment-reservoir-w",          "value":0},     # probability wild reservoir is treated (assume lifelong infection otherwise)
        {"name":"mu-h",             "value":0.000057},  # human mortality rate
        {"name":"mu-r-sd",             "value":0.001369863},    # reservoir mortality rate
        {"name":"mu-r-rd",             "value":0.000176757},    # reservoir mortality rate
        {"name":"mu-r-sw",             "value":0.000156556},    # reservoir mortality rate
        {"name":"mu-r-rw",             "value":0.000182648},    # reservoir mortality rate
        {"name":"wane_immunity",             "value":1/50}    # same for animals and humans
    ]

    def _add_parameter(name: str, value: float):
        model.add_parameter(name, value)

    for parameter in parameters:
        _add_parameter(**parameter)

    # Convenience functions:
    model.add_function("human-population",           "(+ human-susceptible human-exposed human-infectious-one human-infectious-two human-recovered)")
    model.add_function("reservoir-sd-population",       "(+ reservoir-sd-susceptible reservoir-sd-exposed reservoir-sd-infectious reservoir-sd-recovered)")
    model.add_function("reservoir-sw-population",       "(+ reservoir-sw-susceptible reservoir-sw-exposed reservoir-sw-infectious reservoir-sw-recovered)")
    model.add_function("reservoir-rd-population",       "(+ reservoir-rd-susceptible reservoir-rd-exposed reservoir-rd-infectious reservoir-rd-recovered)")
    model.add_function("reservoir-rw-population",       "(+ reservoir-rw-susceptible reservoir-rw-exposed reservoir-rw-infectious reservoir-rw-recovered)")
    model.add_function("tsetse-population",          "(+ tsetse-susceptible tsetse-exposed tsetse-infectious tsetse-non-susceptible)")
    model.add_function("tsetse-human-ratio",         "(/ tsetse-population human-population)")
    model.add_function("tsetse-reservoir-ratio-sd",     "(/ tsetse-population reservoir-sd-population)")
    model.add_function("tsetse-reservoir-ratio-rd",     "(/ tsetse-population reservoir-rd-population)")
    model.add_function("tsetse-reservoir-ratio-sw",     "(/ tsetse-population reservoir-sw-population)")
    model.add_function("tsetse-reservoir-ratio-rw",     "(/ tsetse-population reservoir-rw-population)")
    model.add_function("reservoir-population",       "(+ reservoir-sd-population reservoir-rd-population reservoir-sw-population reservoir-rw-population)")
    model.add_function("reservoir-infectious",       "(+ reservoir-sd-infectious reservoir-rd-infectious reservoir-sw-infectious reservoir-rw-infectious)")
    model.add_function("p-feed-sd",                  "(* (- 1 p-human-feed) (/ reservoir-sd-population (+ reservoir-population non-reservoir-hosts)))")
    model.add_function("p-feed-sw",                  "(* (- 1 p-human-feed) (/ reservoir-sw-population (+ reservoir-population non-reservoir-hosts)))")
    model.add_function("p-feed-rd",                  "(* (- 1 p-human-feed) (/ reservoir-rd-population (+ reservoir-population non-reservoir-hosts)))")
    model.add_function("p-feed-rw",                  "(* (- 1 p-human-feed) (/ reservoir-rw-population (+ reservoir-population non-reservoir-hosts)))")
    model.add_function("p-reservoir-feed",           "(+ p-feed-sd p-feed-sw p-feed-rd p-feed-rw)")

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
    model.add_function("infectious-feed", "(* beta-v (+ (* p-human-feed (/ human-infectious-one human-population)) (* p-reservoir-feed (/ reservoir-infectious reservoir-population))))")#probablility a given feed is infectious and transmission occurs
    model.add_reaction("feed-and-infected",      ["tsetse-susceptible"], ["tsetse-exposed"], "(* p-feed infectious-feed tsetse-susceptible)") #probabilty a feed happens in first 24 hours, is infectious, and transmission occurs
    model.add_reaction("become-non-susceptible", ["tsetse-susceptible"], ["tsetse-non-susceptible"], "(+ (- 1 p-feed)(* p-feed (- 1 infectious-feed)))")#probability tsetse doesn't feed in first 24 hours, or does but transmission doesn't occur
    model.add_reaction("tsetse-progress-to-infectious", ["tsetse-exposed"], ["tsetse-infectious"], "(* sigma-v tsetse-exposed)")

    # domestic swine reservoir S->E->I->R
    model.add_reaction("reservoir-infection-sd",          ["reservoir-sd-susceptible"], ["reservoir-sd-exposed"],    "(* reservoir-sd-susceptible beta-r p-feed p-feed-sd tsetse-reservoir-ratio-sd (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious-sd", ["reservoir-sd-exposed"],     ["reservoir-sd-infectious"], "(* phi-r-s reservoir-sd-exposed)")
    model.add_reaction("reservoir-recovery-sd",           ["reservoir-sd-infectious"],  ["reservoir-sd-recovered"],  "(+ (* omega-r-t treatment-reservoir-d reservoir-sd-infectious) (* omega-r-nt-s (- 1 treatment-reservoir-d) reservoir-sd-infectious))")
    model.add_reaction("reservoir-waning-immunity-sd",    ["reservoir-sd-recovered"],   ["reservoir-sd-susceptible"], "(* wane_immunity reservoir-sd-recovered)")

    #Domestic ruminant reservoir
    model.add_reaction("reservoir-infection-rd",          ["reservoir-rd-susceptible"], ["reservoir-rd-exposed"],    "(* reservoir-rd-susceptible beta-r p-feed p-feed-rd tsetse-reservoir-ratio-rd (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious-rd", ["reservoir-rd-exposed"],     ["reservoir-rd-infectious"], "(* phi-r-r reservoir-rd-exposed)")
    model.add_reaction("reservoir-recovery-rd",           ["reservoir-rd-infectious"],  ["reservoir-rd-recovered"],  "(+ (* omega-r-t treatment-reservoir-d reservoir-rd-infectious) (* omega-r-nt-r (- 1 treatment-reservoir-d) reservoir-rd-infectious))")
    model.add_reaction("reservoir-waning-immunity-rd",    ["reservoir-rd-recovered"],   ["reservoir-rd-susceptible"], "(* wane_immunity reservoir-rd-recovered)")

    #Wild boar reservoir
    model.add_reaction("reservoir-infection-sw",          ["reservoir-sw-susceptible"], ["reservoir-sw-exposed"],    "(* reservoir-sw-susceptible beta-r p-feed p-feed-sw tsetse-reservoir-ratio-sw (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious-sw", ["reservoir-sw-exposed"],     ["reservoir-sw-infectious"], "(* phi-r-s reservoir-sw-exposed)")
    model.add_reaction("reservoir-recovery-sw",           ["reservoir-sw-infectious"],  ["reservoir-sw-recovered"],  "(+ (* omega-r-t treatment-reservoir-w reservoir-sw-infectious) (* omega-r-nt-s (- 1 treatment-reservoir-w) reservoir-sw-infectious))")
    model.add_reaction("reservoir-waning-immunity-sw",    ["reservoir-sw-recovered"],   ["reservoir-sw-susceptible"], "(* wane_immunity reservoir-sw-recovered)")

    #Wild ruminant reservoir
    model.add_reaction("reservoir-infection-rw",          ["reservoir-rw-susceptible"], ["reservoir-rw-exposed"],    "(* reservoir-rw-susceptible beta-r p-feed p-feed-rw tsetse-reservoir-ratio-rw (/ tsetse-infectious tsetse-population))")
    model.add_reaction("reservoir-exposed-infectious-rw", ["reservoir-rw-exposed"],     ["reservoir-rw-infectious"], "(* phi-r-r reservoir-rw-exposed)")
    model.add_reaction("reservoir-recovery-rw",           ["reservoir-rw-infectious"],  ["reservoir-rw-recovered"],  "(+ (* omega-r-t treatment-reservoir-w reservoir-rw-infectious) (* omega-r-nt-r (- 1 treatment-reservoir-w) reservoir-rw-infectious))")
    model.add_reaction("reservoir-waning-immunity-rw",    ["reservoir-rw-recovered"],   ["reservoir-rw-susceptible"], "(* wane_immunity reservoir-rw-recovered)")

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
    # model.add_reaction("reservoir-susceptible-death-birth", ["reservoir-susceptible"], ["reservoir-susceptible"], "(* mu-r reservoir-susceptible)")
    model.add_reaction("reservoir-exposed-death-birth-sd",     ["reservoir-sd-exposed"],     ["reservoir-sd-susceptible"], "(* mu-r-sd reservoir-sd-exposed)")
    model.add_reaction("reservoir-infectious-death-birth-sd",  ["reservoir-sd-infectious"],  ["reservoir-sd-susceptible"], "(* mu-r-sd reservoir-sd-infectious)")
    model.add_reaction("reservoir-recovered-death-birth-sd",   ["reservoir-sd-recovered"],   ["reservoir-sd-susceptible"], "(* mu-r-sd reservoir-sd-recovered)")

    model.add_reaction("reservoir-exposed-death-birth-rd",     ["reservoir-rd-exposed"],     ["reservoir-rd-susceptible"], "(* mu-r-rd reservoir-rd-exposed)")
    model.add_reaction("reservoir-infectious-death-birth-rd",  ["reservoir-rd-infectious"],  ["reservoir-rd-susceptible"], "(* mu-r-rd reservoir-rd-infectious)")
    model.add_reaction("reservoir-recovered-death-birth-rd",   ["reservoir-rd-recovered"],   ["reservoir-rd-susceptible"], "(* mu-r-rd reservoir-rd-recovered)")

    model.add_reaction("reservoir-exposed-death-birth-sw",     ["reservoir-sw-exposed"],     ["reservoir-sw-susceptible"], "(* mu-r-sw reservoir-sw-exposed)")
    model.add_reaction("reservoir-infectious-death-birth-sw",  ["reservoir-sw-infectious"],  ["reservoir-sw-susceptible"], "(* mu-r-sw reservoir-sw-infectious)")
    model.add_reaction("reservoir-recovered-death-birth-sw",   ["reservoir-sw-recovered"],   ["reservoir-sw-susceptible"], "(* mu-r-sw reservoir-sw-recovered)")

    model.add_reaction("reservoir-exposed-death-birth-rw",     ["reservoir-rw-exposed"],     ["reservoir-rw-susceptible"], "(* mu-r-rw reservoir-rw-exposed)")
    model.add_reaction("reservoir-infectious-death-birth-rw",  ["reservoir-rw-infectious"],  ["reservoir-rw-susceptible"], "(* mu-r-rw reservoir-rw-infectious)")
    model.add_reaction("reservoir-recovered-death-birth-rw",   ["reservoir-rw-recovered"],   ["reservoir-rw-susceptible"], "(* mu-r-rw reservoir-rw-recovered)")

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
