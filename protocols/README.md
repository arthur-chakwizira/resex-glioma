# Gradient waveforms: ResEx exchange mapping in glioma

These are the two free gradient waveforms used in *Characterising Water Exchange in Gliomas Using Diffusion MRI with Free Gradient Waveforms* (Chakwizira et al., Imaging Neuroscience). They differ in exchange weighting (Γ ≈ 10 ms vs 40 ms) and have near-identical restriction weighting (V_ω ≈ 3.7·10³ vs 3.6·10³ s⁻²).

| Waveform | Γ [ms] | V_ω [s⁻²] | Peak \|g\| at b_max [mT/m] |
|---|---|---|---|
| 1 | 9.78 | 3746 | 77.8 |
| 2 | 40.16 | 3593 | 72.1 |

## Files

| File | Contents |
|---|---|
| `waveform1.txt`, `waveform2.txt` | Time [s] and effective gradient [T/m] at the maximum b-value (4 ms/µm²). Use these if you only need the waveforms. |
| `glioma_xps.mat` | The full acquisition protocol as a MATLAB struct `xps`, as used for all analyses in the paper (md-dmri conventions). |

## Conventions

- **Effective gradient.** All waveforms are the *effective* gradient g(t): the sign of the gradient after the refocusing pulse is already inverted. They are therefore not the waveforms played out on the scanner, but they can be used directly to compute q(t), b, spectra, etc.
- **Sampling.** 120 samples on a 1 ms raster, covering the 120 ms encoding window. The zero-valued samples near the middle (54–61 ms for waveform 1, 55–64 ms for waveform 2) are the pause in which the refocusing pulse was played. Waveform 1 ends at 115 ms; the remaining samples are zero.
- **Gyromagnetic ratio.** γ = 2π · 42.576·10⁶ rad s⁻¹ T⁻¹. With a slightly different value of γ, recomputed b-values differ from the stored ones at the 10⁻⁵ level.
- **Units.** SI throughout.

## The `xps` struct

Rows are acquisitions (8 in total). The two waveforms alternate, and each was acquired at b = 0, 1.33, 2.67 and 4.0 ms/µm². At each non-zero b-value, 10 directions were acquired and powder-averaged; directions are therefore not included. The waveforms at lower b-values are the b_max waveforms scaled by √(b/b_max), and the b = 0 rows contain all-zero waveforms.

| Field | Size | Unit | Description |
|---|---|---|---|
| `b` | 8×1 | s/m² | b-value |
| `b_delta` | 8×1 | – | b-tensor shape (1 = linear encoding) |
| `gwf` | 8×120 | T/m | effective gradient waveform |
| `dt` | 8×1 | s | time step of `gwf` |
| `Gamma` | 8×1 | s | exchange weighting, Γ = 2∫ t q₄(t) dt |
| `Vomega` | 8×1 | s⁻² | restriction weighting, V_ω = (γ²/b)∫ g(t)² dt |
| `q4` | 8×120 | s⁻¹ | q₄(t) = (1/b²)∫ q²(t′) q²(t′+t) dt′ (zero for b = 0 rows) |
| `s_ind` | 8×1 | – | waveform index (1 or 2) |
| `pa_w` | 8×1 | – | number of directions powder-averaged |
| `n` | 1×1 | – | number of rows |

Here q(t) = γ∫₀ᵗ g(t′) dt′ and b = ∫ q(t)² dt. All integrals are evaluated on the 1 ms raster.

V_ω is the second moment of the two-sided encoding power spectrum, V_ω = (1/2πb)∫ |q(ω)|² ω² dω, which equals the time-domain expression above.

## Quick start

**MATLAB**

```matlab
load('glioma_xps.mat', 'xps')
gamma = 2*pi*42.576e6;
g  = xps.gwf(end-1, :);                  % waveform 1 at b_max
dt = xps.dt(1);
q  = gamma * cumsum(g) * dt;
b  = sum(q.^2) * dt;                     % ~4e9 s/m^2
```

**Python**

```python
import numpy as np
t, g = np.loadtxt('waveform1.txt', unpack=True)
dt = t[1] - t[0]
gamma = 2 * np.pi * 42.576e6
q = gamma * np.cumsum(g) * dt
b = np.sum(q**2) * dt                    # ~4e9 s/m^2
Vomega = gamma**2 * np.sum(g**2) * dt / b

# encoding power spectrum; pad generously, the spectrum is narrow
N = 64 * len(q)
Q2 = np.abs(np.fft.fft(q, N) * dt)**2
f = np.fft.fftfreq(N, dt)                # Hz (omega = 2*pi*f)
pos = f >= 0
centroid_Hz = np.sum(Q2[pos] * f[pos]) / np.sum(Q2[pos])   # ~7.9 Hz (wf 1), ~6.8 Hz (wf 2)
```

## Citation

If you use these waveforms, please cite the paper above.
