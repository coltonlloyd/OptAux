from __future__ import print_function, absolute_import

import numpy as np
import pandas as pd

from cobra.flux_analysis import flux_variability_analysis
import cobra
import cobra.test


def knockout_reactions(model_cons, reaction_list):
    model = model_cons.copy()

    for r in reaction_list:
        model.reactions.get_by_id(r).knock_out()

    return model


def make_consumption_envelope(model, reaction_list, target, num_points=10):
    #model.reactions.get_by_id(target).lower_bound = 0.

    model_ko = knockout_reactions(model, reaction_list)
    model_wt = model.copy()

    max_gr_ko = model_ko.optimize().f
    max_gr_wt = model_wt.optimize().f

    frac_max_list = np.linspace(0, 1, num_points)
    print(len(frac_max_list))

    plot_wt_target = np.zeros(num_points)
    plot_wt_gr = np.zeros(num_points)
    plot_ko_target = np.zeros(num_points)
    plot_ko_gr = np.zeros(num_points)
    index = 0
    for i, frac in enumerate(frac_max_list):
        index += 1 if i % 2 == 0 and i != 0 else 0

        fva_wt = flux_variability_analysis(model_wt, reaction_list=[target],
                                           fraction_of_optimum=frac)
        fva_ko = flux_variability_analysis(model_ko, reaction_list=[target],
                                           fraction_of_optimum=frac)

        plot_wt_target[index] = fva_wt[target]['minimum']
        plot_wt_target[-index-1] = fva_wt[target]['maximum']
        plot_wt_gr[index] = frac * max_gr_wt
        plot_wt_gr[-index-1] = frac * max_gr_wt

        plot_ko_target[index] = fva_ko[target]['minimum']
        plot_ko_target[-index-1] = fva_ko[target]['maximum']
        plot_ko_gr[index] = frac * max_gr_ko
        plot_ko_gr[-index-1] = frac * max_gr_ko

    return ([plot_wt_gr, plot_wt_target]), ([plot_ko_gr, plot_ko_target])
