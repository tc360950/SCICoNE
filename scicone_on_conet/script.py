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


parameters = {
  "bp_detection":
  {
    "window_size":20,
    "threshold":3,
    "bp_limit":200,
  },
  "inference":
  {
    "cluster_trees":
    {
      "n_reps":10,
      "n_iters":4000,
      "move_probs":[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.01, 0.1, 0.01, 1.0, 0.01],
      "n_tries":3
    },
    "full_trees":
    {
      "n_reps":10,
      "n_iters":10000,
      "move_probs":[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.01, 0.1, 0.01, 1.0, 0.01],
    },
  }
}

# RUN SCICoNE
seed = 42  # for reproducibility


sci = scicone.SCICoNE(install_path, temporary_outpath, persistence=False)

# Detect breakpoints
bps = sci.detect_breakpoints(conet_cc, window_size=parameters['bp_detection']['window_size'], threshold=parameters['bp_detection']['threshold'], bp_limit=parameters['bp_detection']['bp_limit'])

print("Breakpoints detected")
scicone.plotting.plot_matrix(conet_cc, bps=bps['segmented_regions'], cbar_title='Raw counts', vmax=4)

# Run cluster trees
alpha = 1./bps['segmented_regions'].shape[0]
gamma = 1./np.mean(np.sum(conet_cc, axis=1))/conet_cc.shape[0] # 1/coverage

sci.learn_tree(conet_cc, bps['segmented_region_sizes'], verbosity=2, n_reps=parameters['inference']['cluster_trees']['n_reps'], cluster=True, full=False, cluster_tree_n_iters=parameters['inference']['cluster_trees']['n_iters'], max_tries=parameters['inference']['cluster_trees']['n_tries'], robustness_thr=0.5, alpha=alpha, max_scoring=True, gamma=gamma)

print(sci.best_full_tree.outputs)
inferred_cns = sci.best_full_tree.outputs['inferred_cnvs']
print(f"RMSE: {sqrt(np.sum((inferred_cns - counts)**2) / counts.size)}")
