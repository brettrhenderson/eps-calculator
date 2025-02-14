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
parser.add_argument('--extend', '-e', action='store_true', help='keep appending steps rather than overwriting them for step 9 (if reset_counters used by accident).')
parser.add_argument('--nosort', '-n', action='store_true', help='do not sort filenames. Use the input order.')
parser.add_argument('--efield', type=float, default=0.001, help='Electric field in au, default=0.001.')
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
if not args.nosort:
    files.sort()

if not args.extend:
    pols = {}
    for file in files:
        print(file)
        _, pol = parsers.get_dipole(file)
        if len(pol):
            for i, p in enumerate(pol):
                pols[pol[i, 0]] = pol[i]
    p2 = OrderedDict(sorted(pols.items()))
    p2 = np.concatenate([pols[step].reshape(1,-1) for step in pols], axis=0)
else:
    p2 = None
    for file in files:
        print(file)
        _, pol = parsers.get_dipole(file)
        if len(pol):
            if p2 is None:
                p2 = pol
            else:
                pol[:, 0] += p2[-1, 0]
                p2 = np.concatenate([p2, pol], axis=0)

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

eps_r = 1 + 4*np.pi * (atot[-1]) / vol / args.efield               # atot already is referenced to initial pol
eps_inf = 1 + 4*np.pi * (p1[-1, -1] - p0[-1]) / vol / args.efield  # assumes no jumps in clamped ion step

print(f'pq = {pquant}, p_i = {p0[-1]}, p_ci = {p1[-1, -1]}, p_f = {p2[-1, -1]}, vol={vol}, efield = {args.efield}')
print(
f"""
Dielectric Constants:

|  High Frequency   |  Low Frequency   | 
|-------------------|------------------|
|        {eps_inf:.3f}      |       {eps_r:.3f}     | 
"""
)

if args.plot:
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
