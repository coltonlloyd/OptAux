from os import system, makedirs
from os.path import join, abspath, dirname


TORQUE_TEMPLATE = """#!/bin/bash
#PBS -q regular
#PBS -l mem=8gb
#PBS -d /home/colton/OptAux/OptAux/submit_jobs
#PBS -e /home/colton/OptAux/OptAux/submit_jobs/error
#PBS -o /home/colton/OptAux/OptAux/submit_jobs/output
"""

CMD_TEMPLATE = """
source /home/colton/.virtualenvs/optaux/bin/activate
python /home/colton/OptAux/OptAux/ME_community/simulate_model.py %s %s %s %s
"""
here = dirname(abspath(__file__))


def apply_torque(pair, f, q1, q2):
    scrip_dir = abspath(join(here, 'commands'))
    job_name = '_'.join([pair, str(f), str(q1), str(q2)])
    job_file = join(scrip_dir, '%s.sh' % job_name)
    with open(job_file, 'w') as outfile:
        outfile.write(TORQUE_TEMPLATE)
        outfile.write(CMD_TEMPLATE % (pair, f, q1, q1))
        outfile.write('\n')
    system('chmod +x ' + job_file)
    system('qsub ' + job_file)
