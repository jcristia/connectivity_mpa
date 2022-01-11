# Do some basic matrix multiplication to figure out what I need to do

import numpy as np

# the proportion of the population that leaves
prop_leave = 0.1


patch_pop = np.array([
    100,
    200,
    300
])


pop_l = patch_pop * prop_leave


disp_prob = np.array([
    [0.1, 0.2, 0.3],
    [0.2, 0.3, 0.1],
    [0.3, 0.4, 0.5]
])


immigrate = np.matmul(disp_prob, pop_l)

# a row is what comes into a patch
# a column is what leaves

# Look at it manually
# patch 1 (10 starting)
# coming in:
1
4
9
# going away:
1
2
3

# It looks like with matrix mult, I get just the total amount coming into each patch.
# so then, sum the columns to get the total proportion that should be leaving
dp_sum = np.sum(disp_prob, axis=0)
emigrate = pop_l * dp_sum
# Yes, this is how I would get the total amount going to other patches,
# BUT, I don't need this number. Where things go from and to gets accounted for
# when I calculate immigration. Then to get the net amount, I just subract
# the total amount that leaves.

# then subtract those and add to the original amount
net = patch_pop - pop_l + immigrate - emigrate
#hmmmm, no. I'm double subtracking with pop_l and emigrate
net = patch_pop - (pop_l - emigrate) + immigrate - emigrate
# but now, my retention is canceling out
retention = pop_l * disp_prob.diagonal()
# so now:
net = patch_pop - (pop_l - emigrate) + immigrate - emigrate + retention

# WRONG
# I need to do this:
net = patch_pop - pop_l + immigrate



####### Make new script from here.
# That is all the logic I need to work through for this script.


##############
# EIGENVECTOR/VALUE STUFF

from numpy import linalg as LA
w, v = LA.eig(disp_prob)

# look at my w and compare them to net
# my patch 1 is fairng the best (largest eigenvalue)
# but, since the other two will go to zero soon, it will also
# therefore, with nothing greater than 1 then this is not a persistent
# metapopulation

# WORKING BACKWARDS TO UNDERSTSAND THE MATH:
# https://web.physics.utah.edu/~detar/lessons/python/numpy_eigen/node1.html
left = np.dot(disp_prob, v[:,0])
right = w[0] * v[:,0]
# so from this, there is an eigenvalue that is representing each row (so the amount coming into each patch).
# This is what eigenvector centrality is. Go look at my old notes in evernote.
# The larger it is the more stuff coming in and the more important it is.
# Like PageRank, a page is important if a lot of stuff points to it.
# So I guess...
# if the leading eigenvalue is greater than 1 then there is positive growth for that patch?
# I don't get how that one value is meant to represent positive growth for the entire pop.
# Hmmmmmm, maybe it is more so about... as long as there is 1 patch with an
# eigenvalue greater than 1 then something is persisting to send shit there, which I guess means that
# the metapopulation is persisting because for 1 patch to persist then at least
# 1 other patch must be persisting?

# http://ordination.okstate.edu/eigen.htm
# Eigenvalues are usually ranked from the greatest to the least. The first eigenvalue is often called the "dominant" or "leading" eigenvalue. Eigenvalues are also often called "latent values".
# The eigenvalue is a measure of the strength of an axis, the amount of variation along an axis, and ideally the importance of an ecological gradient. The precise meaning depends on the ordination method used.

