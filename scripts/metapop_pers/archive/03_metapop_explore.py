# Use simple quantities to test stochastic model

# Based on Patricks metacommunity model:
# https://github.com/plthompson/mcomsimr/blob/master/R/MC_simulate.R

import numpy as np
import pandas as pd


timesteps = 100

# the proportion of the population that leaves at each time step
prop_leave = 0.1

# r naught
r0 = 1.5

# starting populations
patch_pop = np.array([
    100,
    200,
    300
])

# Dispersal probabilities
# Think of it as: a row is the amount that comes into one patch from other
# patches. A column is the amount that leaves a patch to go to the other patches.
disp_prob = np.array([
    [0.1, 0.2, 0.3],
    [0.2, 0.3, 0.1],
    [0.3, 0.4, 0.5]
])

# initialize
pop_net = patch_pop


for i in list(range(timesteps)):

    ### apply logistic growth with Beverton-Holt model ####

    # use starting population as carrying capacity
    # with Beverton-Holt, K = (r0 - 1) * (1/alpha)
    # r0 must be greater than 1
    alpha = (r0 - 1) / patch_pop
    pop_net = pop_net * r0 / (1 + pop_net * alpha)
    pop_net[pop_net<0] = 0

    # so what is Patrick doing in his script?
    # I think it is a lotka voltera competition style thing.
    # Left of the divide sign is the amount that it could grow to. Then on the right
    # that is limited by the amount that is reduced by interspecific competition.
    # He adds 1, otherwise he would be dividing by an amount less than 1 and it would
    # just grow.
    # This is why he said he is setting K implicitly. Pops are limited by
    # competition, and in the previous step, r is set by environmental optima, so 
    # the carrying capacity isn't based on patch size or anything.


    #### expected population size ###
    # add stochasticity to account for deaths (perhaps from enviornmental stochasticity)
    # use poisson distribution to ensure I get a whole number of individuals
    pop_adj = np.random.poisson(pop_net)

    ### amount that disperses ###
    # use a binomial distribution to:
    # (1) add stochasticity around the amount that disperses at each step
    # (2) ensure that I get whole number of individuals dispersing
    # binomial = probability of a 1/0 event
    pop_disp = np.random.binomial(n=pop_adj, p=prop_leave)


    ### Where each individual goes ###

    # approach if we want fractional values:
    #immigrate = np.matmul(disp_prob, pop_disp)

    # approach if we want to move whole individuals:
    # (for each column from a pactch to others, decide where each individual goes
    # based on dispersal probabilities to those patches). This requires going
    # through each column.

    # NUMPY approach
    # add a new row so columns add to 1
    if i==0:
        newrow = 1 - np.sum(disp_prob, axis=0)
        disp_prob = np.vstack([disp_prob, newrow])
    # assign each particle based on dispersal probabilities
    def randomChoice(col):
        samp = np.random.choice(a=disp_prob.shape[0], size=pop_disp[col], p=disp_prob[:,col])
        unique, counts = np.unique(samp, return_counts=True)
        samp = np.zeros(disp_prob.shape[0])
        np.put(samp, unique, counts)
        samp = samp[:-1]
        return(samp)
    # for each column, figure out where it's particles go
    samp_all = map(randomChoice, list(range(disp_prob.shape[1])))
    samp_all = list(samp_all)
    immigrate = sum(samp_all)

    # PANDAS approach
    # newrow = 1 - np.sum(disp_prob, axis=0)
    # disp_prob = np.vstack([disp_prob, newrow])
    # def randomChoice(col):
    #     df = pd.DataFrame(disp_prob[:,col], columns=['probs'])
    #     df['patch'] = df.index
    #     samp = df.sample(n=pop_disp[col], weights='probs', replace=True)
    #     samp = samp.patch.value_counts()
    #     return(samp)
    # samp_all = map(randomChoice, list(range(disp_prob.shape[1])))
    # samp_sum = pd.concat(list(samp_all), axis=1).sum(1).sort_index()
    # samp_sum = samp_sum[:-1] # drop last row
    # immigrate = samp_sum.to_numpy()


    # net
    net = pop_adj - pop_disp + immigrate
    pop_net = net
