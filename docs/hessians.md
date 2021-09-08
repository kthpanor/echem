# Hessians

If we work in AO basis:
```{math}
:label: eq:Lagrangian_in_ao
\frac{\partial L}{\partial \xi}&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{p,q}\omega_{pq}\frac{\partial S_{pq}}{\partial \xi}\\
&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\sum_{p,q}\omega_{pq}C_{\mu p}S^\xi_{\mu\nu}C_{\nu q}\\
&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\omega_{\mu\nu}S^\xi_{\mu\nu}\\
```
and take the derivative:
```{math}
:label: eq:second_order_energy_derivative
\frac{\mathrm{d}^2E}{\mathrm{d}\chi\mathrm{d}\xi}&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial L}{\partial \xi}=\frac{\mathrm{d}}{\mathrm{d}\chi}\left(\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\omega_{\mu\nu}S^\xi_{\mu\nu}\right)\\
&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\frac{\mathrm{d}}{\mathrm{d}\chi}\left(\omega_{\mu\nu}S^\xi_{\mu\nu}\right)\\
&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}\frac{\mathrm{d}}{\mathrm{d}\chi}S^\xi_{\mu\nu}\right)\\
&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
```
By plugging in the expression for the partial derivative of the HF energy in AO basis, we get:
```{math}
:label: eq:second_order_explicit
\frac{\mathrm{d}^2E}{\mathrm{d}\chi\mathrm{d}\xi}&=\frac{\mathrm{d}}{\mathrm{d}\chi}\left( \sum_{\mu \nu} P_{\mu \nu} h_{\mu \nu}^{\xi} + \frac12 \sum_{\mu \nu \lambda \sigma} P_{\mu \nu} P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} + \frac{\partial V_{nn}}{\partial \xi}\right) + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
&= \sum_{\mu \nu} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  h_{\mu \nu}^{\xi} + P_{\mu \nu}h_{\mu\nu}^{\xi\chi}  \right) + \frac12 \sum_{\mu \nu \lambda \sigma} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} + P_{\mu \nu}\frac{\mathrm{d}P_{\lambda \sigma}}{\mathrm{d}\chi}   \langle \mu \lambda || \nu \sigma \rangle^{\xi}+P_{\mu \nu}P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi\chi} \right) \\
&+ \frac{\partial^2 V_{nn}}{\partial \xi\partial \chi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
&= \sum_{\mu \nu} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  h_{\mu \nu}^{\xi} + P_{\mu \nu}h_{\mu\nu}^{\xi\chi}  \right) + \frac12 \sum_{\mu\nu\lambda\sigma}P_{\mu \nu}P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi\chi} + \sum_{\mu \nu \lambda \sigma} \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} +\frac{\partial^2 V_{nn}}{\partial \xi\partial \chi}\\
&+ \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)
```
