from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt

tree = Phylo.read(str(snakemake.input), "newick")
tree.ladderize()

fig = plt.figure(figsize=(13,5), dpi=100)
matplotlib.rc("font", size=8)
matplotlib.rc("xtick", labelsize=12)
matplotlib.rc("ytick", labelsize=12)
axes = fig.add_subplot(1,1,1)
Phylo.draw(tree, axes=axes, do_show=False)
fig.savefig(str(snakemake.output))
