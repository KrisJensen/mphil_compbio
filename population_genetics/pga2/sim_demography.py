#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate msprime with population structure
"""

import numpy as np
import matplotlib.pyplot as plt
import msprime
import time
import pickle
import math


def fullsim(tgen = 25, T_EU_AS = 50000):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_AS = 500
    N_EU = 500
    # Times are provided in years, so we convert into generations.

    T_EU_AS /= tgen
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0#0.004
    r_AS = 0#0.0055
    #N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    #N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.

    m_EU_AS = 0
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=34, initial_size=N_EU),
        msprime.PopulationConfiguration(
            sample_size=4, initial_size=N_AS)]
        #msprime.PopulationConfiguration(
        #    sample_size=1, initial_size=N_AS, growth_rate=r_AS)]
        
    migration_matrix = [[0, m_EU_AS],
                        [m_EU_AS, 0]]
    
    demographic_events = [msprime.MassMigration(
            time=T_EU_AS, source=1, destination=0, proportion=1.0)]
    
    '''
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    '''
    
    replicates = msprime.simulate(
        recombination_rate = 0.5e-9 * tgen,
        length = 100000,
        #length = 2.8*10**9,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        num_replicates = 2)
    
    ts = []
    
    for r in replicates:
        tmins = []
        for tree in r.trees():
            tmin = min([tree.time(tree.mrca(37, 1)),
                        tree.time(tree.mrca(36, 1)),
                        tree.time(tree.mrca(35, 1)),
                        tree.time(tree.mrca(34, 1)),
                        ])
            tmins.append(tmin)
        print('new tmin:', min(tmins))
        ts.append(min(tmins))
        
        l = r.aslist()
        for tree in l:
            print(tree.draw(format="unicode"))
    
    
    #dd = msprime.DemographyDebugger(
    #    population_configurations=population_configurations,
    #    migration_matrix=migration_matrix,
    #    demographic_events=demographic_events)
    #dd.print_history()
    
    return replicates, ts

r, ts = fullsim()
#for l in r:
#    print(l.first().draw(format="unicode"))


