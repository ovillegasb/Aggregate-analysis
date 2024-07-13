#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
This program extracts the coordinates of the smallest and largest aggregates by number of molecules
for each selected frame. The file "cluster_resid.dat" and "clusters.dat" generated by aggregate.py
must be in the run folder.

"""

import pandas as pd
import argparse
import mdtraj as md
import numpy as np
import time


def get_options():
    """ get options from the command line """

    parser = argparse.ArgumentParser(
        prog="Frames",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s [-options]",
        epilog="Enjoy the program!",
        description=__doc__
    )

    parser.add_argument("-f",
                        help="Trajectory file.",
                        metavar="FILE",
                        default="traj_comp.xtc",
                        type=str)

    parser.add_argument("-s",
                        help=("Topology file."),
                        metavar="GRO",
                        default="top.gro",
                        type=str)

    parser.add_argument("-n", "--nmol",
                        help=("Number of molecules in the system."),
                        metavar="NMOL",
                        default=5,
                        type=int)

    parser.add_argument("-dt", "--timestep",
                        help=("Number of step by time unit."),
                        metavar="TIMESTEP",
                        default=20,
                        type=int)

    parser.add_argument("-out",
                        help=("Output directory."),
                        metavar="OUTPUT",
                        default=".",
                        type=str)

    return vars(parser.parse_args())


def center_of_mass(coords, masses, L):
    """ compute the center of mass, the mass weighterd barycenter """
    theta_i = coords / L * 2 * np.pi
    xi_i = np.cos(theta_i)
    eta_i = np.sin(theta_i)
    xi_m = np.sum(xi_i * masses[:, np.newaxis], axis=0) / masses.sum()
    eta_m = np.sum(eta_i * masses[:, np.newaxis], axis=0) / masses.sum()
    theta_m = np.arctan2(-eta_m, -xi_m) + np.pi

    # return self.wcoords.sum(axis=0) / self.masses.sum()
    return L * theta_m / 2 / np.pi


def cluster_centering(coords, masses, L):
    """ completes the broken cluster molecules due to the pbc condition. """

    # Traslation using center of mass like reference poit
    newcoords = list()
    cm = center_of_mass(coords, masses, L)
    for i in coords:
        r_i = list()
        for j in range(len(i)):
            if np.abs(i[j] - cm[j]) < L / 2:
                r_i.append(i[j])
            elif np.abs(i[j] - L - cm[j]) < L / 2:
                r_i.append(i[j] - L)
            elif np.abs(i[j] + L - cm[j]) < L / 2:
                r_i.append(i[j] + L)
        newcoords.append(np.array(r_i))

    return np.array(newcoords)


def extract_frames(file, nmol, dt):
    """ Analysis aggregates """
    header = ['Rg1', 'Rg2', 'Rg3', 'Shape', 'v', 'mm', 'd', 'rho', 'nc', 'frame']
    clust = pd.read_csv(
        file,
        sep=r'\s+',
        skiprows=13,
        names=header
    )
    clust['frame'] *= 1
    clust['time'] = clust['frame'] * dt

    # ifram   nmol   [resid]
    resid = []
    with open('cluster_resid.dat') as RESID:
        for line in RESID:
            if '#' not in line:
                resid.append([str(i) for i in line.split()[2:]])

    clust["resid"] = resid

    frame_clust = []

    for n in range(1, nmol + 1):
        a = clust[clust.nc == n]
        a = a[a.d == a.d.max()]
        frame_clust.append(a)

    for n in range(1, nmol + 1):
        a = clust[clust.nc == n]
        a = a[a.d == a.d.min()]
        frame_clust.append(a)

    frames = pd.concat(frame_clust, ignore_index=True)
    frames = frames.sort_values(by=['frame'], ignore_index=True)
    return frames


def save_clusters_coords(frames, traj, out="."):

    for i in frames.index:
        resid = frames.loc[i, 'resid']
        resid = ' '.join(resid)
        ndx = traj[frames.loc[i, 'frame']].top.select(f'resid {resid}')
        aggr = traj[frames.loc[i, 'frame']].atom_slice(ndx)

        # aggr.save_gro(f"{out}/conf{frames.loc[i, 'frame']}_{frames.loc[i, 'nc']}_{i}.gro")

        # gmx trjconv -s run.tpr -f confout.gro -pbc cluster -o newcoords.gro

        coords = aggr.xyz[0]
        masses = np.array([atom.element.mass for atom in aggr.topology.atoms])
        L = aggr.unitcell_lengths[0][0]

        newcoords = cluster_centering(coords, masses, L)

        aggr.xyz = np.array([newcoords])

        # saving to gro file.
        aggr.save_gro(f"{out}/aggr{frames.loc[i, 'frame']}_{frames.loc[i, 'nc']}_{i}.gro")

    lines = ''
    for i in frames.index:
        lines += '%s    ' % frames.loc[i, 'frame']
        resid = frames.loc[i, 'resid']
        lines += '-'.join(resid)
        lines += '\n'

    with open(f'{out}/frames.dat', 'w') as f:
        f.write(lines)


def main():
    # starting setup
    args = get_options()
    t0 = time.time()

    nmol = args['nmol']
    dt = args['timestep']

    traj = args['f']
    top = args['s']

    # output directory
    out = args['out']

    print('\nSearching frames.')
    print('%s molecules' % nmol)
    frames = extract_frames('clusters.dat', nmol, dt)
    traj = md.load(traj, top=top)

    print("Saving coordinates.")
    save_clusters_coords(frames, traj, out)
    dt = time.time() - t0
    print("Done in %.0f s\n" % dt)


if __name__ == '__main__':
    main()
