# smpesgm2017
Smart meter data analysis for the conference paper submitted to PESGM2018

This is the algorithm and experimentation code accompanying the paper:

Yanchao Liu, Tingli Hu, Caisheng Wang. Timely Detection of Abnormal Inactivity Using Smart Meter Data. 
http://yliu.eng.wayne.edu/research/occupancy.pdf

## Abstract
A novel algorithm is proposed in this paper for detecting
abnormal inactivities within a single-occupied household
based on the smart meter readings. Such inactivities include
immobilizing medical conditions or sudden deaths of elderly
or disabled occupants who live alone, the delayed discovery of
which poses realistic social concerns as the population ages. By
extracting novel features from the power variation and employing
the state-of-the-art probabilistic methods for anomaly detection,
the algorithm is able to cold start from limited historical data and
perform well without extended parameter tuning. The method
is validated on a real data set with simulated scenarios and is
shown to be effective. The algorithm is small profile in data usage
and therefore has potential to be implemented in a distributed
and embedded system, such as in a smart meter itself.

## Algorithm Highlights
1. Defined predictive features that reflect diurnal patterns and emerging trends in power consumption within an individual household.

2. Applied Grubb's test based on (novel) weighted Mahalanobis distance.

3. Requires little training data to start working
