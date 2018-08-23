from __future__ import absolute_import, division, print_function
from optaux.me_community.me_model_community import make_binary_community_me
from optaux import resources
import pickle
from os.path import abspath, dirname

here = dirname(abspath(__file__))
resource_dir = resources.__path__[0]

with open('%s/iJL1678b.pickle' % resource_dir, 'rb') as f:
    me = pickle.load(f)

with open('%s/iJL1678b.pickle' % resource_dir, 'rb') as f:
    me_cons = pickle.load(f)

me.unmodeled_protein_fraction = .8
me_cons.unmodeled_protein_fraction = .8

print("USING 65 for ALL keffs")
for r in me.reactions:
    if hasattr(r, 'keff'):
        r.keff = 65.
for r in me.process_data:
    if hasattr(r, 'keff'):
        r.keff = 65.

for r in me_cons.reactions:
    if hasattr(r, 'keff'):
        r.keff = 65.
for r in me_cons.process_data:
    if hasattr(r, 'keff'):
        r.keff = 65.


me.global_info['k_deg'] = 0.
me.update()
me_cons.global_info['k_deg'] = 0.
me_cons.update()
print('remove rna degradation')

make_binary_community_me(me, me_cons)
