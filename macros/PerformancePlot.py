from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import statistics

squrt_E = [0.687,0.687,0.687,0.783,0.529,0.529,0.529,0.529,0.529,0.707,0.707,0.707,0.613,0.613,0.612,0.612]
sigma_E_over_E = [2.8005,2.8590,2.8392,3.9356,1.5165,1.5106,1.5042,1.5000,1.5352,3.1485,3.0055,3.1350,2.1613,2.2000,2.1944,2.4258]
uncert_sigma_E_over_E = [0.01132,0.01792,0.00896,0.01227,0.00363,0.00754,0.00307,0.00307,0.00307,0.05350,0.07900,0.02750,0.03684,0.04023,0.10337,0.02996]

plt.errorbar(squrt_E, sigma_E_over_E, xerr=None, yerr=uncert_sigma_E_over_E, fmt='o', capsize=5)

## Polynomial fit (degree=2)
#coeffs = np.polyfit(squrt_E, sigma_E_over_E , 2)
#poly = np.poly1d(coeffs)
#x_fit = np.linspace(0.400, 0.950, 16)
#y_fit = poly(x_fit)

## Quad Sum Fit
def quad_poly_func(x, c0, c1, c2):
    return np.sqrt( np.square(c0 * np.square(x)) + np.square(c1 * x) + np.square(c2) )

popt, pcov = curve_fit(quad_poly_func, squrt_E, sigma_E_over_E)
c0_fit, c1_fit, c2_fit = popt

x_fit = np.linspace(0.400, 0.950, 16)
y_fit = quad_poly_func(x_fit, c0_fit, c1_fit, c2_fit)

## Define chi2/ndf
#chi_sum = 0
#n = 0
#for i in squrt_E:
#    fit_value = coeffs[0] * np.square(i) + coeffs[1] * i + coeffs[2]
#    chi_sqr = np.square(sigma_E_over_E[n]-fit_value) / fit_value
#    chi_sum += chi_sqr
#    n+=1
#ndf = len(squrt_E) - (2 - 1)
#goodFit = chi_sum / ndf

## Define R^2
resid_sum = 0
mean_sum = 0
mean_y = statistics.mean(sigma_E_over_E)
n = 0
for i in squrt_E:
    fit_value = np.sqrt( np.square(popt[0] * np.square(i)) + np.square(popt[1] * i) + np.square(popt[2]) )
    
    resid_sq = np.square(sigma_E_over_E[n] - fit_value)
    resid_sum += resid_sq

    mean_sq = np.square(sigma_E_over_E[n] - mean_y)
    mean_sum += mean_sq

    n+=1

r_sq = 1 - (resid_sum / mean_sum)

## Plotting
plt.plot(x_fit, y_fit, color='green', label=f'Quadrature Sum Fit\ny = ' + "{:.3f}".format(popt[0]) + 'x$^2$ \u2295 ' + "{:.7f}".format(popt[1]) + 'x \u2295 ' + "{:.8f}".format(popt[2]) + "\nR$^2$ = " +"{:.3f}".format(r_sq))
plt.legend()

plt.xlabel("1/\u221A(E) [GeV$^{-1/2}$]")
plt.ylabel("\u03C3$_E$/E [%]")
plt.title("BBCal Energy Resolution Performance")

plt.show()
