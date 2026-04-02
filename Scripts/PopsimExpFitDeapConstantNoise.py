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
    parser.add_argument("--EXPR_MEAN", type=float, required=True,
                        help="Mean expression of initial population")
    parser.add_argument("--EXPR_SD", type=float, required=True,
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
# Main simulation
# =====================================================
def main(args):

    # -----------------------------
    # Reproducibility
    # -----------------------------
    random.seed(42)
    np.random.seed(42)

    EXPR_MIN = 0.01
    EXPR_MAX = 2.0

    # -----------------------------
    # Fitness functions
    # -----------------------------
    if args.fitness_function == "gaussian":
        pdf_grid = norm.pdf(
            np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01),
            loc=args.FIT_var1,
            scale=args.FIT_var2,
        )
        c = 112 / pdf_grid.max()
        #c = 56 / pdf_grid.max()
        #c = 90 / pdf_grid.max()

        def doubling_time(x):
            #return 240 - c * norm.pdf(
            #    x, loc=args.FIT_var1, scale=args.FIT_var2
            #)
            return 191 - c * norm.pdf(
                x, loc=args.FIT_var1, scale=args.FIT_var2
            )
            
    elif args.fitness_function == "lognorm":
        pdf_grid = lognorm.pdf(
            np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01),
            scale=np.exp(args.FIT_var1), 
            loc=0, 
            s= args.FIT_var2
        )
        c = 90 / pdf_grid.max()

        def doubling_time(x):
            return 160 - c * lognorm.pdf(
                x, scale=np.exp(args.FIT_var1), loc=0, s= args.FIT_var2
            )

    elif args.fitness_function == "mixednorm":
        weights = [args.weight, (1-args.weight)]
        EXP = np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01)
        pdf_grid = (weights[0] * norm.pdf(EXP, loc=args.FIT_var1, scale=args.FIT_var2) + 
             weights[1] * norm.pdf(EXP, loc=args.FIT_var1_pair, scale=args.FIT_var2_pair))/2

        c = 90 / pdf_grid.max()


        def doubling_time(x):
                    return 160 - c * ((weights[0] * norm.pdf(x, loc=args.FIT_var1, scale=args.FIT_var2) + 
                     weights[1] * norm.pdf(x, loc=args.FIT_var1_pair, scale=args.FIT_var2_pair))/2)


    elif args.fitness_function == "mixedlognorm":
        weights = [args.weight, (1-args.weight)]
        EXP = np.arange(EXPR_MIN, EXPR_MAX + 0.01, 0.01)
        pdf_grid = (weights[0] * lognorm.pdf(EXP, scale=np.exp(args.FIT_var1), s=args.FIT_var2, loc=0) + 
             weights[1] * lognorm.pdf(EXP, scale=np.exp(args.FIT_var1_pair), s=args.FIT_var2_pair, loc=0))/2

        c = 90 / pdf_grid.max()


        def doubling_time(x):
                    return 160 - c * ((weights[0] * lognorm.pdf(x, scale=np.exp(args.FIT_var1), s=args.FIT_var2, loc=0) + 
                     weights[1] * lognorm.pdf(x, scale=np.exp(args.FIT_var1_pair), s=args.FIT_var2_pair, loc=0))/2)


        
    else:
        raise ValueError(f"Unknown fitness function: {args.fitness_function}")

    # -----------------------------
    # DEAP setup
    # -----------------------------
    if not hasattr(creator, "FitnessMax"):
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register(
        "expression", random.gauss, args.EXPR_MEAN, args.EXPR_SD
    )
    toolbox.register(
        "individual",
        tools.initRepeat,
        creator.Individual,
        toolbox.expression,
        n=1
    )
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("select", tools.selRoulette)

    def fitness(ind):
        return (1.0 / doubling_time(ind[0]),)

    toolbox.register("evaluate", fitness)

    # -----------------------------
    # Run simulations
    # -----------------------------
    records = []

    for it in range(1, args.iterations + 1):
        population = toolbox.population(n=args.pop_size)

        for ind in population:
            ind.fitness.values = toolbox.evaluate(ind)

        time = 0.0

        while time <= args.total_time:
            expressions = np.array([ind[0] for ind in population])

            mean_expr = expressions.mean()
            var_expr = expressions.var()
            sd_expr = np.sqrt(var_expr)

            dt_vals = doubling_time(expressions)
            mean_dt = dt_vals.mean()

            cv_expr = sd_expr / mean_expr if mean_expr > 0 else np.nan
            fano_expr = var_expr / mean_expr if mean_expr > 0 else np.nan

            records.append({
                "Initial_ExpMean": args.EXPR_MEAN,
                "Initial_ExpSD": args.EXPR_SD,
                "Fitfun": args.fitness_function,
                "FIT_var1": args.FIT_var1,
                "FIT_var2": args.FIT_var2,
                "Heritability": args.h,
                "iteration": it,
                "time": time,
                "mean_doubling_time": mean_dt,
                "SD": sd_expr,
                "Fano": fano_expr,
                "CV": cv_expr,
            })

            # Selection
            mothers = toolbox.select(population, k=args.pop_size)

            # Reproduction with heritability
            offspring = []
            for mom in mothers:
                mom_expr = mom[0]

                daughter_expr = (
                    args.EXPR_MEAN
                    + np.sqrt(args.h) * (mom_expr - args.EXPR_MEAN) #changed to the defined initial mean, not "mean_expr"
                    + np.sqrt(1 - args.h) * np.random.normal(0, args.EXPR_SD) #changed to the defined initial noise, not "sd_expr"
                )

                daughter_expr = np.clip(
                    daughter_expr, EXPR_MIN, EXPR_MAX
                )

                child = creator.Individual([daughter_expr])
                child.fitness.values = toolbox.evaluate(child)
                offspring.append(child)

            population = offspring
            time += args.dt

    # -----------------------------
    # Save CSV
    # -----------------------------
    os.makedirs(args.output_dir, exist_ok=True)

    #filename = (
    #    f"SIM_mean_{args.EXPR_MEAN}_"
    #    f"sd_{args.EXPR_SD}_"
    #    f"her_{args.h}_"
    #    f"{args.fitness_function}.csv"
    #)

    filepath = os.path.join(args.output_dir, args.outfile)
    pd.DataFrame.from_records(records).to_csv(filepath, index=False)

    print(f"\nSimulation complete.")
    print(f"Results saved to:\n{filepath}")


# =====================================================
# Entry point
# =====================================================
if __name__ == "__main__":
    args = parse_args()
    main(args)
