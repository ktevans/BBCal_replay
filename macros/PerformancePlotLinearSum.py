from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import statistics

## Define dataset from results
energy = [2.12,1.63,3.58,2.00,2.66,2.67] #using weighted averages ## E, scattered electron energy as defined by BB tracking/optics, treated here as a constant value with no uncertainty
squrt_E = [0.68680,0.78326,0.52852,0.70711,0.61314,0.61199] #using weighted averages ## 1/sqrt(E) 
sigma_E_over_E = [2.83698,3.93558,1.51280,3.07557,2.18149,2.24640] #using weighted averages ## sigma_E is the sigma on the fit of E/p for different calibration results
uncert_sigma_E_over_E = [0.01765,0.01227,0.00538,0.04904,0.01934,0.09661] #using weighted averages ## the uncertainty on sigma_E is found from the E/p fit made during calibrations

## Make scatterplot with error bars
plt.errorbar(energy, sigma_E_over_E, xerr=None, yerr=uncert_sigma_E_over_E, fmt='.', capsize=4)

## Define a polynomial fit where y is (sigma/E) and x is (E)
def poly_func(x, c0, c1, c2):
    return (c0 / np.sqrt(x)) + (c1 / x) + c2 ## sigma_E/E = (c0/sqrt(E)) + (c1/E) + c2

## Create a curve fit based on the defined poly_fit
popt, pcov = curve_fit(poly_func, energy, sigma_E_over_E, sigma=uncert_sigma_E_over_E) ## note that curve_fit uses a least squares fitting algorithm
c0_fit, c1_fit, c2_fit = popt
#print("\ncovariance matrix:\n")
#print(pcov)

c0_uncert = (np.sqrt(pcov[[0],[0]]))
c1_uncert = np.sqrt(pcov[[1],[1]])
c2_uncert = np.sqrt(pcov[[2],[2]])

print("Fit parameter values:")
print("{:.3f}".format(popt[0]) + " \u00B1 " + "{:.3f}".format(c0_uncert[0]) + "%")
print("{:.3f}".format(popt[1]) + " \u00B1 " + "{:.3f}".format(c1_uncert[0]) + "%")
print("{:.4f}".format(popt[2]) + " \u00B1 " + "{:.4f}".format(c2_uncert[0]) + "%")


## Define axis and fit limits based on dataset
min_e_lim = np.min(energy) - (0.5 * np.std(energy))
max_e_lim = np.max(energy) + np.std(energy)

## Apply limits to the range of the fit and then use then define the fit for our data
x_fit = np.linspace(min_e_lim, max_e_lim, len(energy)*10)
y_fit = poly_func(x_fit, c0_fit, c1_fit, c2_fit)

## Define R^2
resid_sum = 0 ## necessary variable to create sums
mean_sum = 0 ## necessary variable to create sums
mean_y = statistics.mean(sigma_E_over_E) ## necessary variable to create R^2
n = 0 ## index tracker
for i in energy:
    fit_value = (popt[0] / np.sqrt(i)) + (popt[1] / i) + popt[2] # find y value based on the fit, y_fit
    
    resid_sq = np.square(sigma_E_over_E[n] - fit_value) # resid^2 = (y_true - y_fit)
    resid_sum += resid_sq # add to sum of resid^2

    mean_sq = np.square(sigma_E_over_E[n] - mean_y) # mean^2 = (y_true - mean_y)
    mean_sum += mean_sq # add to sum of mean^2

    n+=1 # move to next index in lists

r_sq = 1 - (resid_sum / mean_sum) # R^2 = 1 - [ sum((y_true-y_fit)^2) / sum((y_true - mean_y)^2) ]

## Plotting data and fit with fit equation and R^2 in legend
plt.plot(x_fit, y_fit, color='green', label=f'Second Order Polynomial Fit\n\u03C3$_E$/E\'$_e$ = ' + "{:.3f}".format(popt[0]) + '% /\u221A(E\') + ' + "{:.4f}".format(popt[2]) + '% + ' + "{:.3f}".format(popt[1]) + '%/E\'' + "\nR$^2$ = " +"{:.3f}".format(r_sq))
plt.legend()
plt.xlabel("E'$_e$ [GeV]")
plt.ylabel("\u03C3$_{E'}$/E'$_e$ [%]")
plt.title("BBCal Energy Resolution Performance")

plt.show()
