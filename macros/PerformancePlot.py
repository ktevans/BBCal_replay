from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import statistics

squrt_E = [0.68680,0.78326,0.52852,0.70711,0.61314,0.61199] #using weighted averages
#[0.687,0.687,0.687,0.783,0.529,0.529,0.529,0.529,0.529,0.707,0.707,0.707,0.613,0.613,0.612,0.612]

energy = [2.12,1.63,3.58,2.00,2.66,2.67] #using weighted averages
#[2.120,2.120,2.120,1.630,3.580,3.580,3.580,3.580,3.580,2.000,2.000,2.000,2.660,2.660,2.670,2.670]

sigma_E_over_E = [2.83698,3.93558,1.51280,3.07557,2.18149,2.24640] #using weighted averages
#[2.8005,2.8590,2.8392,3.9356,1.5165,1.5106,1.5042,1.5000,1.5352,3.1485,3.0055,3.1350,2.1613,2.2000,2.1944,2.4258]

uncert_sigma_E_over_E = [0.01765,0.01227,0.00538,0.04904,0.01934,0.09661] #using weighted averages
#[0.01132,0.01792,0.00896,0.01227,0.00363,0.00754,0.00307,0.00307,0.00307,0.05350,0.07900,0.02750,0.03684,0.04023,0.10337,0.02996]

#plt.errorbar(squrt_E, sigma_E_over_E, xerr=None, yerr=uncert_sigma_E_over_E, fmt='o', capsize=5)
plt.errorbar(energy, sigma_E_over_E, xerr=None, yerr=uncert_sigma_E_over_E, fmt='.', capsize=4)

## Polynomial fit (degree=2)
#coeffs = np.polyfit(squrt_E, sigma_E_over_E , 2, w=uncert_sigma_E_over_E)
#poly = np.poly1d(coeffs)
#x_fit = np.linspace(1.20, 3.60, 16)
#y_fit = poly(x_fit)

## poly Fit
#def poly_func(x, c0, c1, c2):
def poly_func(x, c0, c2):
    return np.sqrt( np.square(c0 * (1/np.sqrt(x))) + np.square(c2) )
    #return np.sqrt( np.square(c0 * (1/np.sqrt(x))) + np.square(c1 * (1/x)) + np.square(c2) )
    #return (c0 / np.sqrt(x)) + (c1 / x) + c2

popt, pcov = curve_fit(poly_func, energy, sigma_E_over_E, sigma=uncert_sigma_E_over_E)
#c0_fit, c1_fit, c2_fit = popt
c0_fit, c2_fit = popt
print(pcov)

min_e_lim = np.min(energy) - (0.5 * np.std(energy))
max_e_lim = np.max(energy) + np.std(energy)

x_fit = np.linspace(min_e_lim, max_e_lim, len(energy)*10)
#y_fit = poly_func(x_fit, c0_fit, c1_fit, c2_fit)
y_fit = poly_func(x_fit, c0_fit, c2_fit)

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
for i in energy:
    fit_value = np.sqrt( np.square(popt[0] * (1/np.sqrt(i))) + np.square(popt[1]) )
    #fit_value = np.sqrt( np.square(popt[0] * (1/np.sqrt(i))) + np.square(popt[1] * (1/i)) + np.square(popt[2]) )
    #fit_value = coeffs[0] * np.square(i) + coeffs[1] * i + coeffs[2]
    #fit_value = (popt[0] / np.sqrt(i)) + (popt[1] / i) + popt[2]
    
    resid_sq = np.square(sigma_E_over_E[n] - fit_value)
    resid_sum += resid_sq

    mean_sq = np.square(sigma_E_over_E[n] - mean_y)
    mean_sum += mean_sq

    n+=1

r_sq = 1 - (resid_sum / mean_sum)

## Plotting
#plt.plot(x_fit, y_fit, color='green', label=f'Second Order Polynomial Fit\n\u03C3$_E$/E\'$_e$ = ' + "{:.3f}".format(popt[0]) + '% /\u221A(E\') + ' + "{:.4f}".format(popt[2]) + '% + ' + "{:.3f}".format(popt[1]) + '%/E\'' + "\nR$^2$ = " +"{:.3f}".format(r_sq))
#plt.plot(x_fit, y_fit, color='green', label=f'Quadrature Sum Fit\n\u03C3$_E$/E\'$_e$ = ' + "{:.6f}".format(popt[0]) + '% /\u221A(E\') \u2295 ' + "{:.8f}".format(popt[2]) + '% \n\u2295 ' + "{:.3f}".format(popt[1]) + '%/E\'' + "\nR$^2$ = " +"{:.3f}".format(r_sq))
plt.plot(x_fit, y_fit, color='green', label=f'Quadrature Sum Fit\n\u03C3$_E$/E\'$_e$ = ' + "{:.6f}".format(popt[0]) + '% /\u221A(E\') \u2295 ' + "{:.8f}".format(popt[1]) + "%\nR$^2$ = " +"{:.3f}".format(r_sq))
plt.legend()
#\u2295


#plt.xlabel("1/\u221A(E) [GeV$^{-1/2}$]")
plt.xlabel("E'$_e$ [GeV]")
plt.ylabel("\u03C3$_{E'}$/E'$_e$ [%]")
plt.title("BBCal Energy Resolution Performance")

plt.show()
