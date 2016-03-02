# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:43:57 2015

@author: lic14
"""

import numpy as np
from scipy.stats import rv_discrete as rvd
from scipy.stats import gaussian_kde

def weight(log_weight):
    '''convert log_likelihood to weight'''
    nsam = np.shape(log_weight)[0]
    weight = np.zeros((nsam,))
    for j in range(nsam):
        weight[j] = 1.0/np.sum(np.exp(log_weight-log_weight[j]))
    return weight

def PF(weights, prior):
    '''weights: (n,) array; prior: (n, d) array; where n is the number of sample and d is number of nodes'''
    (n,_) = np.shape(prior)
    post_indices = rvd(name='post_indices',values = (np.arange(n),weights)) # Resample particles according to their weights # generate a random variables named post_indices
    post_samples = post_indices.rvs(size=n) # generate random numbers for this variables.

    pos = prior[post_samples]
    return pos
 
def resample(sample):
    '''sample: 1-d array. Resample the featrue by ksdensity.
    Note: plese only apply this to continuous variables'''

	# resample by kernel density estimate
    kde = gaussian_kde(sample)
    MAX = np.max(sample)
    MIN = np.min(sample)
    (nsam, ) = np.shape(sample)
    new_sample = np.array([])
    while np.shape(new_sample)[0] < nsam - 2:
        new_sample_sub = kde.resample()
        new_sample_sub = np.ravel(new_sample_sub)
        new_sample_sub = new_sample_sub[new_sample_sub >= MIN]
        new_sample_sub = new_sample_sub[new_sample_sub <= MAX]
        new_sample = np.concatenate((new_sample, new_sample_sub))
        
    new_sample = new_sample[0:nsam-2]
    new_sample = np.concatenate((new_sample, np.array([MIN, MAX])))
        

    # rearrange of the new samples at the same rank as the old sample, so that the correlation with other varialbes are maintained
    sort_index = np.argsort(sample)
    p_d1_sample_sort = np.sort(new_sample)
    new_sample_ordered = np.zeros(len(new_sample))
    j = 0
    for i in sort_index:
        new_sample_ordered[i] = p_d1_sample_sort[j]
        j = j+1
    return new_sample_ordered