import numpy as np
import random
import scicone

from math import sqrt

from conet.generative_model import CNSampler, EventTree, CountsGenerator, EventTreeGenerator

# SCICoNE path
install_path = '/code/build/' # must be a full path
temporary_outpath = './tmp'


loci = 1500
tree_size = 10
cells = 200
random.seed(2243534)
np.random.seed(233454)

# generate event tree and cell data, this might take a while
cn_s = CNSampler.create_default_sampler()
t_gen = EventTreeGenerator(cn_sampler=cn_s, tree_size=tree_size, no_loci=loci)
tree: EventTree = t_gen.generate_random_tree()
d_gen = CountsGenerator(cn_s, tree)
counts, attachment, corrected_counts, brkp_matrix = d_gen.generate_data(loci, cells)
conet_cc = np.transpose(np.array(corrected_counts)[:, 5:])



# RUN SCICoNE
seed = 42  # for reproducibility


sci = scicone.SCICoNE(install_path, temporary_outpath, persistence=False)

# Detect breakpoints
bps = sci.detect_breakpoints(conet_cc, threshold=3.0)
print("Breakpoints detected")
scicone.plotting.plot_matrix(conet_cc, bps=bps['segmented_regions'], cbar_title='Raw counts', vmax=4)

sci.learn_tree(conet_cc, bps['segmented_region_sizes'], n_reps=4, seed=seed, cluster=True, full=True)

print(sci.best_full_tree.outputs)
inferred_cns = sci.best_full_tree.outputs['inferred_cnvs']
print(f"RMSE: {sqrt(np.sum((inferred_cns - counts)**2) / counts.size)}")
