from __future__ import print_function, division, absolute_import

from multiprocessing import Pool
import os

from optaux import me_community


sim_script_dir = me_community.__path__[0]


def run_pool(function, values, processes=4):
    pool = Pool(processes=processes)
    pool.map(function, values)


def submit_job(value_tuple):
    pair, f, q1, q2, mode, restrict, k1, k2, uptake, model = value_tuple
    os.system("python3.6 %s/simulate_model.py %s %s %s %s %s "
              "--Restrict_crossfeeding %s --keff_transporter_1 %s "
              "--keff_transporter_2 %s --glucose_uptake %s "
              "--model_to_use %s" %
              (sim_script_dir, pair, f, q1, q2, mode, restrict, k1, k2,
               uptake, model))


def does_sim_exist(simulation_directory, pair, mode, f, q1, q2, k1, k2,
                   met=None):

    if mode == 'secretion_keff_sweep':
        file_loc = \
            '%s/%s/%s/%.2f_%.2f_0_secretion_multiplier/' \
            '%.2f_frac_strain1_sol.json' % \
            (simulation_directory, mode, pair, k1*100, k2*100, f)

    elif mode == 'unmodeled_sweep':
        file_loc = \
            '%s/%s/%s/%.2f_%.2f_0_unmodeled_protein/' \
            '%.2f_frac_strain1_sol.json' % \
            (simulation_directory, mode, pair, q1*100, q2*100, f)

    elif mode == 'metabolite_limitation':
        file_loc = \
            '%s/%s/%s/%s/%.2f_frac_strain1_sol.json' % \
            (simulation_directory, mode, pair, met, f)

    elif mode == 'default':
        file_loc = '%s/%s/%s/%.2f_%.2f_unmodeled_protein' \
                   '/%.2f_frac_strain1_sol.json' % \
                   (simulation_directory, mode, pair, q1*100, q2*100, f)
    elif mode == 'glucose_limited':
        file_loc = '%s/%s/%s/%.2f_%.2f_unmodeled_protein' \
                   '/%.2f_frac_strain1_sol.json' % \
                   (simulation_directory, mode, pair, q1*100, q2*100, f)
    else:
        raise UserWarning('Mode (%s) is not valid' % mode)

    if os.path.exists(file_loc):
        return True
    else:
        print('Output in ', file_loc)
        return False
