# Special aspects

## Comparison to experiment

- Beer-Lambert
- Excitation: vertical, 0-0, adiabatic
- Franck-Condon, Herzberg-Teller, ...
- Selection rules

Here the overlap of the S<sub>0</sub> vibrational ground state and the third vibrational state of S<sub>i</sub> is such that this transition will contribute more to the intensity than, for instance, vibrational ground state to vibrational ground state. Calculating the different contributions to the different ground states (the Franck–Condon factors) thus yields a good first improvement of the theoretical spectra, featuring a smoother spectrum. If the PES of S<sub>0</sub> and S<sub>i</sub> are sufficiently close, the inclusion of only ground state to ground state may be sufficient, designated 0–0 in the figure. Experimentally, it may only be possible to resolve the 0–0 transition, or the wavelength corresponding to the maximum of absorption, $\lambda$<sub>max</sub> . The vertical transition energy thus gives only a rough estimate of the actual physical situation, but it is in most cases sufficiently accurate.

## Vibrations

## Broadening

The choice of convolution function depends on the dominant broadening effect (lifetime of excited/ionized state, experimental resolution, Doppler effect, ...). We generally focus on Lorentzian broadening, which is the type which naturally comes out of damped response calculations and which relates to the finite life-time of the excited/ionized state.

Giving each ionization energy equal weight, Lorentzian and Gaussian broadening are performed as below, for which the full-width at half-max (FWHM) is equal to $\gamma $ and $ 2\sqrt{2\ln{2}} \sigma$, respectively. If better correspondence to experiment is sought, a Voigt profile can be constructed from the convolution of a Lorentzian and a Gaussian.

```python
def lorentzian(x, y, xmin, xmax, xstep, gamma):
    '''
    Lorentzian broadening function
    
    Call: xi,yi = lorentzian(energies, intensities, start energy, end energy, energy step, gamma)
    '''
    xi = np.arange(xmin,xmax,xstep); yi=np.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] = yi[i] + y[k] * gamma / ( (xi[i]-x[k])**2 + (gamma/2.)**2 ) / np.pi
    return xi,yi

def gaussian(x, y, xmin, xmax, xstep, sigma):
    '''
    Gaussian broadening function
    
    Call: xi,yi = gaussian(energies, intensities, start energy, end energy, energy step, gamma)
    '''
    xi = np.arange(xmin,xmax,xstep); yi=np.zeros(len(xi))
    for i in range(len(xi)): 
        for k in range(len(y)): yi[i] = yi[i] + y[k]*np.e**(-((xi[i]-x[k])**2)/(2*sigma**2))
    return xi,yi

# IEP:s from above calculation
esca_ies = [302.187, 299.362, 295.405, 293.277]

plt.figure(figsize=(10,5))
plt.subplot(221); plt.title('Lorentzian, FWHM = 0.4 eV')
x,y = esca_ies,np.ones((len(esca_ies)))
xi,yi = lorentzian(x,y,min(x)-2,max(x)+2,0.01,0.4); plt.plot(xi,yi)

plt.subplot(222); plt.title('Gaussian, FWHM = 0.4 eV')
xi,yi = gaussian(x,y,min(x)-2,max(x)+2,0.01,0.4/(2*np.sqrt(2*np.log(2)))); plt.plot(xi,yi)

plt.subplot(223); plt.title('Both, peak max normalized')
xi,yi = gaussian(x,y,min(x)-2,max(x)+2,0.01,0.4/(2*np.sqrt(2*np.log(2)))); plt.plot(xi,yi/max(yi))
xi,yi = lorentzian(x,y,min(x)-2,max(x)+2,0.01,0.4); plt.plot(xi,yi/max(yi))

plt.subplot(224); plt.title('Both, area normalized')
xi,yi = gaussian(x,y,min(x)-2,max(x)+2,0.01,0.4/(2*np.sqrt(2*np.log(2)))); plt.plot(xi,yi/sum(yi))
xi,yi = lorentzian(x,y,min(x)-2,max(x)+2,0.01,0.4); plt.plot(xi,yi/sum(yi))
plt.tight_layout(); plt.show()
```

<!-- #region -->


## Charge-transfer
<!-- #endregion -->
