# Optimized-Confidence-Interval
This simple program would help you to get the confidence interval of some asymmetric distributions, such as F distribution and Chi-Square distribution.

# Introduction
During interval estimation, there would be varies of distributions pivot variable, some have excellent features, such as normal distrubution, but others may not be that awesome.
Chi-Square distribution and F distribution are two distributions that the pivot variable usually have, however, they're not symmetric, which means if you take the upper a / 2 and uppser 1 - a / 2 quantile of the distribution as the interval of your pivot variable, it would not be the interval with minimum length, though it's convenient.
The program can generate the approximately optimized interval of the pivot variable, given the distribution and the confidence level of the pivot variable.

# Algorithm
The method of obtaining the interval is based on composite Simpson's rule for integration and used given values for the calculation of Gamma function. For more details, you can refer to the Thesis in the folder.
