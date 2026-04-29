#!/usr/bin/env python3

import os
import argparse
import random
import numpy as np
import pandas as pd
from scipy.stats import norm, lognorm
from deap import base, creator, tools


# =====================================================
# Argument parsing
# =====================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate expression-based evolution with heritability"
    )
    parser.add_argument("--EXPR_MEAN_A", type=float, required=True,
                        help="Mean expression of initial population")
    parser.add_argument("--EXPR_SD_A", type=float, required=True,
                        help="Expression SD of initial population")
    parser.add_argument("--EXPR_MEAN_B", type=float, required=True,
                        help="Mean expression of initial population")
    parser.add_argument("--EXPR_SD_B", type=float, required=True,
                        help="Expression SD of initial population")
    parser.add_argument("--FIT_var1", type=float, required=True,
                        help="Mean (optimal) expression for fitness function")
    parser.add_argument("--FIT_var1_pair", type=float, required=False,default=None,
                        help="Mean (2nd optimal) expression for fitness function")
    parser.add_argument("--FIT_var2", type=float, required=True,
                        help="SD of expression in fitness function")
    parser.add_argument("--FIT_var2_pair", type=float, required=False,default=None,
                        help="SD of expression (2nd optimal) in fitness function")
    parser.add_argument("--weight", type=float, default=1,
                        help="Weight of the first distribution for mixed distribution sims; Max value=1, Min value =0")
    parser.add_argument("--h", type=float, required=True,
                        help="Heritability (0–1)")
    parser.add_argument("--iterations", type=int, default=20,
                        help="Number of replicate simulations")
    parser.add_argument("--pop_size", type=int, default=10000,
                        help="Population size")
    parser.add_argument("--total_time", type=float, default=1000.0,
                        help="Total simulation time (minutes)")
    parser.add_argument("--dt", type=float, default=20.0,
                        help="Time step (minutes)")
    parser.add_argument("--fitness_function", type=str, default="gaussian",
                        help="Fitness function type (currently: gaussian, lognorm, mixednorm, mixedlognorm )")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory for output CSV")
    parser.add_argument("--outfile", type=str, required=True,
                        help="Name for output CSV")

    return parser.parse_args()


# =====================================================
# Fitness functions
# =====================================================
def build_doubling_time(args, EXPR_MIN, EXPR_MAX):

    if args.fitness_function == "gaussian":
        pdf_grid = norm.pdf(
            np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01),
            loc=args.FIT_var1,
            scale=args.FIT_var2,
        )
        c = 112 / pdf_grid.max()

        def doubling_time(x):
            return 191 - c * norm.pdf(x, loc=args.FIT_var1, scale=args.FIT_var2)

    elif args.fitness_function == "lognorm":
        pdf_grid = lognorm.pdf(
            np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01),
            scale=np.exp(args.FIT_var1),
            loc=0,
            s=args.FIT_var2,
        )
        c = 90 / pdf_grid.max()

        def doubling_time(x):
            return 160 - c * lognorm.pdf(
                x, scale=np.exp(args.FIT_var1), loc=0, s=args.FIT_var2
            )

    elif args.fitness_function == "mixednorm":
        weights = [args.weight, 1 - args.weight]
        EXP = np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01)

        pdf_grid = (
            weights[0] * norm.pdf(EXP, args.FIT_var1, args.FIT_var2)
            + weights[1] * norm.pdf(EXP, args.FIT_var1_pair, args.FIT_var2_pair)
        ) / 2

        c = 90 / pdf_grid.max()

        def doubling_time(x):
            return 160 - c * (
                weights[0] * norm.pdf(x, args.FIT_var1, args.FIT_var2)
                + weights[1] * norm.pdf(x, args.FIT_var1_pair, args.FIT_var2_pair)
            ) / 2

    elif args.fitness_function == "mixedlognorm":
        weights = [args.weight, 1 - args.weight]
        EXP = np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01)

        pdf_grid = (
            weights[0] * lognorm.pdf(EXP, s=args.FIT_var2, scale=np.exp(args.FIT_var1))
            + weights[1]
            * lognorm.pdf(EXP, s=args.FIT_var2_pair, scale=np.exp(args.FIT_var1_pair))
        ) / 2

        c = 90 / pdf_grid.max()

        def doubling_time(x):
            return 160 - c * (
                weights[0]
                * lognorm.pdf(x, s=args.FIT_var2, scale=np.exp(args.FIT_var1))
                + weights[1]
                * lognorm.pdf(
                    x, s=args.FIT_var2_pair, scale=np.exp(args.FIT_var1_pair)
                )
            ) / 2

    else:
        raise ValueError("Unknown fitness function")

    return doubling_time


# =====================================================
# Main simulation
# =====================================================
def main(args):

    random.seed(42)
    np.random.seed(42)

    EXPR_MIN = 0.01
    EXPR_MAX = 2.0

    doubling_time = build_doubling_time(args, EXPR_MIN, EXPR_MAX)

    # -----------------------------
    # DEAP setup
    # -----------------------------
    if not hasattr(creator, "FitnessMax"):
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register("select", tools.selRoulette)

    def fitness(ind):
        return (1.0 / doubling_time(ind[0]),)

    toolbox.register("evaluate", fitness)

    # -----------------------------
    # Simulation
    # -----------------------------
    records = []

    for it in range(1, args.iterations + 1):

        # Initialize two populations
        population = []

        for _ in range(args.pop_size):
            ind = creator.Individual([random.gauss(args.EXPR_MEAN_A, args.EXPR_SD_A)])
            ind.label = "A"
            population.append(ind)

        for _ in range(args.pop_size):
            ind = creator.Individual([random.gauss(args.EXPR_MEAN_B, args.EXPR_SD_B)])
            ind.label = "B"
            population.append(ind)

        # Evaluate initial fitness
        for ind in population:
            ind.fitness.values = toolbox.evaluate(ind)

        time = 0.0

        while time <= args.total_time:

            expressions = np.array([ind[0] for ind in population])
            labels = np.array([ind.label for ind in population])

            def stats(mask):
                subset = expressions[mask]
                if len(subset) == 0:
                    return np.nan, np.nan
                return subset.mean(), subset.std()
            
            dt_vals = doubling_time(expressions)

            def mean_dt_for(mask):
                subset = dt_vals[mask]
                return subset.mean() if len(subset) > 0 else np.nan

            mean_A, sd_A = stats(labels == "A")
            mean_B, sd_B = stats(labels == "B")

            freq_A = np.mean(labels == "A")
            freq_B = np.mean(labels == "B")

            mean_dt_A = mean_dt_for(labels == "A")
            mean_dt_B = mean_dt_for(labels == "B")

            PopSize = len(population)

            #mean_dt = doubling_time(expressions).mean()

            records.append({
                "Initial_ExpMean": args.EXPR_MEAN_A,
                "Initial_ExpSD": args.EXPR_SD_A,
                "Fitfun": args.fitness_function,
                "FIT_var1": args.FIT_var1,
                "FIT_var2": args.FIT_var2,
                "Heritability": args.h,
                "iteration": it,
                "time": time,
                "mean_A": mean_A,
                "SD_A": sd_A,
                "mean_B": mean_B,
                "SD_B": sd_B,
                "freq_A": freq_A,
                "freq_B": freq_B,
                "PopSize": PopSize,
                #"mean_doubling_time": mean_dt,
                "mean_dt_A": mean_dt_A,
                "mean_dt_B": mean_dt_B,
            })

            # -----------------------------
            # Selection (competition)
            # -----------------------------
            mothers = toolbox.select(population, k=2 * args.pop_size)

            # -----------------------------
            # Reproduction
            # -----------------------------
            offspring = []

            for mom in mothers:
                mom_expr = mom[0]

                # Use population-specific SD
                if mom.label == "A":
                    sd = args.EXPR_SD_A
                    mean = args.EXPR_MEAN_A
                else:  # "B"
                    sd = args.EXPR_SD_B
                    mean = args.EXPR_MEAN_B

                # Each mother produces TWO daughters
                for _ in range(2):
                    daughter_expr = (
                        mean
                        + np.sqrt(args.h) * (mom_expr - mean)
                        + np.sqrt(1 - args.h) * np.random.normal(0, sd)
                    )

                daughter_expr = np.clip(daughter_expr, EXPR_MIN, EXPR_MAX)

                child = creator.Individual([daughter_expr])
                child.label = mom.label
                child.fitness.values = toolbox.evaluate(child)

                offspring.append(child)


            population = offspring
            time += args.dt

    # -----------------------------
    # Save output
    # -----------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    filepath = os.path.join(args.output_dir, args.outfile)

    pd.DataFrame.from_records(records).to_csv(filepath, index=False)

    print("\nSimulation complete.")
    print(f"Results saved to:\n{filepath}")


# =====================================================
# Entry point
# =====================================================
if __name__ == "__main__":
    args = parse_args()
    main(args)
