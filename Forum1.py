from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt

# using the turtle DNA and after using the muscle alignment tool suing fasta format
with open('alignment_sequence.txt', "r") as aln:
    alignment = AlignIO.read(aln, "fasta")
print(type(alignment))

# calculating the genetic distance
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)
print(distance_matrix)

# make the tree
constructor = DistanceTreeConstructor(calculator)
turtle_tree = constructor.build_tree(alignment)
turtle_tree.rooted = True
print(turtle_tree)

# Make a better looking tree using the features of matplotlib
fig = plt.figure(figsize=(13, 5), dpi=100) # create figure & set the size
matplotlib.rc('font', size=12)              # fontsize of the leaf and node labels
matplotlib.rc('xtick', labelsize=10)       # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=10)       # fontsize of the tick labels
turtle_tree.ladderize()
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(turtle_tree, axes=axes)
fig.savefig("turtles_cladogram")

# write the tree in a file
Phylo.write(turtle_tree, "turtle_tree.xml", "phyloxml")

