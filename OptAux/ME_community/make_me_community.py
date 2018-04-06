from __future__ import absolute_import, division, print_function
from OptAux.ME_community.me_model_community import make_binary_community_me
from OptAux import resources
import pickle
from os.path import abspath, dirname

here = dirname(abspath(__file__))
resource_dir = '%s/resources' % resources.__path__[0]

with open('%s/iJL1678b.pickle' % resource_dir, 'rb') as f:
    me = pickle.load(f)

with open('%s/iJL1678b.pickle' % resource_dir, 'rb') as f:
    me_cons = pickle.load(f)

me.unmodeled_protein_fraction = .8
me_cons.unmodeled_protein_fraction = .8

me.global_info['k_deg'] = 0.
me.update()
me_cons.global_info['k_deg'] = 0.
me_cons.update()
print('remove rna degradation')

make_binary_community_me(me, me_cons)
