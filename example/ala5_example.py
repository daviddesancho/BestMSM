
# coding: utf-8

# # A simple working example of the BestMSM package

# The following provides a minimal example to how the BestMSM package works. We construct an MSM using data from explicit water atomistic simulations of the ala5 peptide, the small molecule that you can see below.
# <img src="files/ala5_nowater.png" width="200">.

# ### Loading the trajectory data

# The first thing that you must do is to import the trajectory module from the BestMSM package.

# In[1]:

import bestmsm.trajectory as traj


# The only required input is the name of the file with strings corresponding to time-stamps and state names (or indexes or whatever you want; they'll be handled as strings). In the current example states look like the following:

# ```1 11111
# 2 11111
# 3 11111
# 4 11111
# 5 11111
# 6 11111
# 7 11111
# 8 11110
# 9 11100
# ...
# ```

# These strings are just helix-coil states for the ala5 pentapeptide. Ones ("1") indicate that a residue is in the helical configuration of the Ramachandran map, while zeros ("0") correspond to the coil region (most of what we see here are fully helical states).

# Then we generate an instance of the TimeSeries class

# In[2]:

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")


# The TimeSeries class has a number of atributes, which we can explore. These attributes are the time stamps of the snapshots of the trajectory ('`time`'), the corresponding the states ('`states`'), and the names of the states ('`keys`'), the lag between snapshots ('`dt`') and the filename ('`filename`'). For example, in this case, the lag between snapshots is 1 unit of time (in this case, picoseconds).

# In[3]:

traj_ala5.dt


# ### Generating an MSM

# Having read one or multiple trajectories we now invoke the MasterMSM class.

# In[4]:

import bestmsm.msm as msm
msm_ala5 = msm.MasterMSM([traj_ala5])


# The instance of the `MasterMSM class` starts off only with the trajectories (or instances of the `trajectory` class) we created above.

# In[5]:

print msm_ala5.keys


# For multiple files, MasterMSM will merge the keys from the different trajectories and make a common set. Alternatively, one can provide the list of states to be considered. This option may be useful in different situations. For example, we may want states to be ordered in a predefined way, or we may want to entirely disregard states that we know are very rarely populated. In this case, the instance of the MasterMSM class would be generated in the following way.

# In[6]:


# This way we have given a list of 32 states which are sorted by the number of helical residues. You can check that the states are ordered so that the first few are mostly coil and the last few are fully helical. This may turn out to be more convenient in this case.

# In[9]:



# Now that you have generated your MSM, you can start computing interesting things. You probably want to generate a count matrix, a transition matrix and a rate matrix, and possibly derive relaxation times and equilibrium probabilities from them. In order to do this, we generate an instance of the `MSM` class. The way things are organized in BestMSM, you may never have to deal with the MSM class directly. Everything can be orchestrated from `MasterMSM`.

# So what we do is to generate the MSM for a given lag time using the data loaded before. Here we naively choose 5 picoseconds

# In[10]:

msm_ala5.do_msm(5)


# When invoking the `do_msm` method a lot of stuff has happened under the hood. First, a transition count matrix has been computed, and then the connectivity of that matrix has been checked. Finally, the transition matrix and its eigenvalues and eigenvectors have been produced.

# ### Chapman-Kolmogorov tests

# Although the example above gives some idea of the available functionalitites, you may want to start finding out which lag time is a sensible choice for your system. In order to do that, you want to know the dependence of your model with respect to the lag time used to construct it. In order to identify the right lag time, we usually look at the convergence of certain properties as a function of the lag time.

# One of the first tests is the dependence of the relaxation times with the lag time used to compute the transition count matrix. This test is readily available in the MasterMSM class.

# In[ ]:

msm_ala5.convergence_test(plot=True, N=5, error=False, sliding=True)


# As you can see, there is a considerable amount of output here. In fact all we are doing is to construct the MSM at different values of the lag time, calculating errors and then plotting the dependence of the relaxation times with the lag. 
# 
# The options we have introduced are `plot=True` and `N=5`, for plotting the 5 slowest eigenmodes after the calculation is done, and `sliding=True` for the sliding window method to be used in the estimation of transition count matrices (see documentation). In the plot the different colours correspond to different eigenvalues. They seem pretty well converged at growing lag times, which suggests that the Markov assumption is acceptable. In this case, for simplicity we have not calculated errors but these should be included for properly identifying the right lag time to use.

# The other usual convergence test for MSMs is an explicit comparison between the simulation data and the model predictions (this test is, in fact, usually called the Chapman-Kolmogorov test). In order to do this, one simulates the relaxation of population from a state or set of states and compares the decay of the population with the actual result from the transition count matrix. This test is also available in BestMSM.

# In[ ]:

msm_ala5.chapman_kolmogorov(init=['00000', '00001', '00010', '00100', '10000', '01000'], \
        sliding=True, plot=True)


msm_ala5.chapman_kolmogorov(init=['11111', '11110', '01111'], \
        sliding=True, plot=True)

## Here what we see as circles is the result from the simulation data and as lines we show the results for a given range of lag times. It can be seen  
#
## ## Working with one MSM
#
## After carrying out this test, you usually choose a lag time for which the MSM seems well converged and then move ahead with the analysis. In our case we are going to use the MSM at a lag time of 5. The MSM at different lag times can be accessed from the MasterMSM instance very simply:
#
## In[13]:
#
#msm_ala5.do_msm(5)
#
#
## Again, we recover the result shown before, but now we can access and manipulate the results.
#
## In[15]:
#
#msm_ala5.msms[5]
#
#
## As you can see msm_ala5.msms[5] is an instance of the MSM class. The trick is that a MasterMSM can keep simultaneously multiple MSM instances, corresponding to different lag times.
#
## ## Error analysis 
#
## One nice feature of the implementation is that we can straightforwardly carry out a bootstrap analysis.
#
## In[18]:
#
#tau_ave, tau_std, peq_ave, peq_std = msm_ala5.msms[5].boots()
#
#
## Depending on the number of available processors, this may take some time. The number of bootstrap samples can be defined by the user, as well as whether results are being plotted.  
#
## ## Committors and fluxes
#
## Another usual result we are interested in is the values of the pfold or committor, or the folding rates from one set of microstates (UU) to another set of microstates (FF). In order to do this calculation we use the Berezhkovskii-Hummer-Szabo method (J Chem Phys, 2009). In our case, we do the calculation as:
#
## In[19]:
#
#J, pfold, kf = msm_ala5.msms[5].do_pfold(FF=31, UU=0)
#
#
## ## PCCA clustering
#
## The time-scale separation between modes can be exploited to generate more intuitive models where each individual state is formed by a cluster of mcirostates. To perform this clustering we import another module.
#
## In[22]:
#
##import bestmsm.pcca import PCCA
#
#
## In[ ]:
#
#
#
