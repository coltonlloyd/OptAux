from __future__ import absolute_import, division, print_function
from optaux.me_community.me_model_community import make_binary_community_me
from optaux import resources
import pickle
from os.path import abspath, dirname
from copy import deepcopy

from cobrame.io.json import load_json_me_model

here = dirname(abspath(__file__))
resource_dir = resources.__path__[0]

print('remove rna degradation')


# Make default community model
me = load_json_me_model('%s/iJL1678b.json' % resource_dir)

me.unmodeled_protein_fraction = .8
me.global_info['k_deg'] = 0.
me.update()

me_cons = deepcopy(me)

make_binary_community_me(me, me_cons, 'community_me_default_keffs.pickle')


# Make community model with all keffs = 65
me = load_json_me_model('%s/iJL1678b.json' % resource_dir)

print("USING 65 for ALL keffs")
for r in me.reactions:
    if hasattr(r, 'keff'):
        r.keff = 65.
for d in me.process_data:
    if hasattr(d, 'keff'):
        d.keff = 65.

for r in me_cons.reactions:
    if hasattr(r, 'keff'):
        r.keff = 65.
for d in me_cons.process_data:
    if hasattr(d, 'keff'):
        d.keff = 65.

me.unmodeled_protein_fraction = .8
me.global_info['k_deg'] = 0.
me.update()

me_cons = deepcopy(me)

make_binary_community_me(me, me_cons, 'community_me_65_keffs.pickle')


# Make community model with all null kappmax keffs
me = load_json_me_model('%s/iJL1678b_null_keffs.json' % resource_dir)
me.unmodeled_protein_fraction = .8
me.global_info['k_deg'] = 0.
me.update()

me_cons = deepcopy(me)

make_binary_community_me(me, me_cons, 'community_me_null_keffs.pickle')
