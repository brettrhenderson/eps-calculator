import argparse
from pybec import parsers, analysis, plotters, utils
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Calculate the dielectric permittivity of a bulk material from Quantum\n" +
                                             "Espresso polarization output files.")
parser.add_argument('zero_field', help='The file containing the zero field polarization output.')
parser.add_argument('relaxed_ion', help='The file or files containing the relaxed ion polarization output.', nargs='+')
parser.add_argument('--clamped-ion', '-c', help='The file containing the clamped ion polarization output.')
parser.add_argument('--plot', '-p', action='store_true', help='Plot the extracted polarization curve to check it.')
parser.add_argument('--efield', type=float, default=0.001, help='Electric field in au, default=0.001.')
parser.add_argument('--eps-bulk', type=float, default=9.26, help='Rel Permittivity of the bulk, for calculating enhancement, default=9.26.')
parser.add_argument('--eps-inf-bulk', type=float, default=3.04, help='High frequency Rel Permittivity of the bulk, default=3.04')
parser.add_argument('--inclusion-element', default='Ag', help='Inclusion element symbol, default=Ag')
args = parser.parse_args()

# Load Cell as a numpy array from an output file
cell = parsers.get_final_cell(args.zero_field)

# volume of cell used for calculating polarization density
vol = analysis.volume(cell)

# for orthorhombic cell, polarization quantum is the length of the z principle axis
pquant = cell[-1, -1]

# Get the last step of zero field polarization
p0 = parsers.get_dipole(args.zero_field)[1]
p0 = p0[-1]

# clamped ion polarization (optional)
if args.clamped_ion:
    p1 = parsers.get_dipole(args.clamped_ion)[1]

# Get the relaxed ion polarization from all steps
files = args.relaxed_ion
files.sort()
pols = {}
for file in files:
    _, pol = parsers.get_dipole(file)
    if len(pol):
        for i, p in enumerate(pol):
            pols[pol[i, 0]] = pol[i]
p2 = OrderedDict(sorted(pols.items()))
p2 = np.concatenate([pols[step].reshape(1,-1) for step in pols], axis=0)

# sum up cumulatively over timesteps and make sure end of step 8 aligns with beginning of 9
if args.clamped_ion:
    # stitch the times from clamped and relaxed ion together
    t8 = np.concatenate([np.array([0]), np.cumsum(p1[:, 1])])  # include last timestep from step 7
    t9 = np.cumsum(p2[:, 1]) + p2[0, 0] * p2[0, 1] + t8[-1]
    time = np.concatenate([t8, t9])

    # Stitch the polarizations together
    atot_uncorrected = np.concatenate([np.array([p0[-1]]), p1[:, -1], p2[:, -1]]) - p0[-1]
    ae_uncorrected = np.concatenate([np.array([p0[2]]), p1[:, 2], p2[:, 2]]) - p0[2]
    ai_uncorrected = np.concatenate([np.array([p0[3]]), p1[:, 3], p2[:, 3]]) - p0[3]

else:
    time = np.concatenate([np.array([0]), np.cumsum(p2[:, 1])])
    atot_uncorrected = np.concatenate([np.array([p0[-1]]), p2[:, -1]]) - p0[-1]
    ae_uncorrected = np.concatenate([np.array([p0[2]]), p2[:, 2]]) - p0[2]
    ai_uncorrected = np.concatenate([np.array([p0[3]]), p2[:, 3]]) - p0[3]

# correct jumps by adding polarization quantum when needed
atot = analysis.correct_jumps(atot_uncorrected, jump_quantum=pquant)
ae = analysis.correct_jumps(ae_uncorrected, jump_quantum=pquant, jump_threshold=10)
ai = analysis.correct_jumps(ai_uncorrected, jump_quantum=pquant)

coords = parsers.get_final_positions(args.zero_field)
al = len(coords[args.inclusion_element])/(sum([len(coords[el]) for el in coords]))
eps_r = 1 + 4*np.pi * (atot[-1]) / vol / args.efield               # atot already is referenced to initial pol
eps_inf = 1 + 4*np.pi * (p1[-1, -1] - p0[-1]) / vol / args.efield  # assumes no jumps in clamped ion step
alpha_r = (eps_r - args.eps_bulk) / (4 * np.pi * al)
alpha_inf = (eps_inf - args.eps_inf_bulk) / (4 * np.pi * al)

print(
f"""
Dielectric Constants:

|  High Frequency  |  Low Frequency  |  alpha_inf  |   alpha_r  | 
|------------------|-----------------|-------------|------------|
|        {eps_inf:.2f}      |       {eps_r:.2f}     |    {alpha_inf:.2f}     |     {alpha_r:.2f}   | 

* alpha values calculated using a bulk matrix relative permittivity of {args.eps_bulk} 
  and high frequency permittivity {args.eps_inf_bulk}.
"""
)

# plot to check that polarization quantum jumps have been smoothed
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,6))

# uncorrected
ax1.set_title('Uncorrected Cell Dipole', fontsize=18)
ax1.plot(time, atot_uncorrected, label="Total")
ax1.plot(time, ae_uncorrected, label="Electronic")
ax1.plot(time, ai_uncorrected, label="Ionic")
ax1.legend(fontsize=16)
ax1.set_xlabel('Time / au', fontsize=16)
ax1.set_ylabel('Cell Dipole', fontsize=16)
ax1.hlines(0, 0,6400, 'k', linestyle='--')

# corrected
ax2.set_title('Corrected Cell Dipole', fontsize=18)
ax2.plot(time, atot, label='Total')
ax2.plot(time, ae, label='Electronic')
ax2.plot(time, ai, label='Ionic')
ax2.hlines(0, 0,6400, 'k', linestyle='--')
ax2.legend()

plt.tight_layout()
plt.show()
