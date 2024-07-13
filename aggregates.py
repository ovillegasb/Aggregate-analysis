#!/usr/bin/env python
# -*- coding=utf-8 -*-

r"""

             _____   _____  _____   ______  _____         _______  ______   _____
     /\     / ____| / ____||  __ \ |  ____|/ ____|    /\ |__   __||  ____| / ____|
    /  \   | |  __ | |  __ | |__) || |__  | |  __    /  \   | |   | |__   | (___
   / /\ \  | | |_ || | |_ ||  _  / |  __| | | |_ |  / /\ \  | |   |  __|   \___ \
  / ____ \ | |__| || |__| || | \ \ | |____| |__| | / ____ \ | |   | |____  ____) |
 /_/    \_\\______| \_____||_|  \_\|______|\_____|/_/    \_\|_|   |______||_____/


Program that performs different analyses to aggregates contained in a system. It takes as main
input files the file to detect the topology a .gro and as trajectory .xtc, .trr or .gro.


"""

import argparse
import time
import numpy as np
from scipy.constants import N_A as avogadro
from scipy.spatial.distance import cdist
import mdtraj as md
import itertools as it

TITLE = r"""

             _____   _____  _____   ______  _____         _______  ______   _____
     /\     / ____| / ____||  __ \ |  ____|/ ____|    /\ |__   __||  ____| / ____|
    /  \   | |  __ | |  __ | |__) || |__  | |  __    /  \   | |   | |__   | (___
   / /\ \  | | |_ || | |_ ||  _  / |  __| | | |_ |  / /\ \  | |   |  __|   \___ \
  / ____ \ | |__| || |__| || | \ \ | |____| |__| | / ____ \ | |   | |____  ____) |
 /_/    \_\\______| \_____||_|  \_\|______|\_____|/_/    \_\|_|   |______||_____/


"""


class GyrationTensor:
    """ Class to represent a Gyration Tensor and compute related quantities """

    def __init__(self, coords, masses, box):
        """
        Args:
            matrix (3x3): 3x3 matrix of the gyration tensor
        """

        """

        new = list()
        for xyz in coords:

            xyz = np.array(xyz)

            new.append(xyz)

        coords = np.array(new)

        """
        self.coords = coords

        self.masses = masses
        # Assuming cubic box
        self.L = box[0][0]

        # Traslation using center of mass like reference poit

        newcoords = list()
        center_of_mass = self.center_of_mass

        for i in coords:
            r_i = list()
            for j in range(len(i)):
                if np.abs(i[j] - center_of_mass[j]) < self.L / 2:
                    r_i.append(i[j])
                elif np.abs(i[j] - self.L - center_of_mass[j]) < self.L / 2:
                    r_i.append(i[j] - self.L)
                elif np.abs(i[j] + self.L - center_of_mass[j]) < self.L / 2:
                    r_i.append(i[j] + self.L)
            newcoords.append(r_i)

        newcoords = np.array(newcoords)
        self.coords = newcoords

        # weighted coords
        self.wcoords = self.coords * self.masses[:, np.newaxis]

    @property
    def center_of_mass(self):
        """ compute the center of mass, the mass weighterd barycenter """

        theta_i = self.coords / self.L * 2 * np.pi

        """
            coords = self.coords.copy()

            print(coords)
            print(type(coords))

            # correct = lambda x: np.array(x)

            # print(np.apply_along_axis(np.array, 0, coords))

            new = list()
            for xyz in coords:
                # print(xyz)
                # print(type(xyz))

                xyz = np.array(xyz)

                # print(xyz)
                # print(type(xyz))

                new.append(xyz)

            coords = np.array(new)

            theta_i = coords / self.L * 2 * np.pi

            print(theta_i)

            exit()
            # for i in range(len(self.coords)):
            #     self.coords[i] = np.array(self.coords[i], dtype=np.float64)

            # theta_i = self.coords / self.L * 2 * np.pi

        """

        """
            print('TypeError: Check types')
            print(self.coords)
            print(type(self.coords))
            print(self.coords.dtype)

            # print(self.coords.astype(np.ndarray))
            # print(self.coords.astype(np.ndarray).dtype)

            print(self.coords[0])
            print(type(self.coords[0]))
            print(np.array(self.coords[0]).dtype)

            print(self.coords.astype(np.ndarray))

            coords = np.array(self.coords, dtype=np.ndarray)
            print(coords)
            exit()
        """
        # else:
        #     print(TypeError)
        #     print("Problems calculing center og mass")
        #     exit()

        xi_i = np.cos(theta_i)
        eta_i = np.sin(theta_i)
        xi_m = np.sum(xi_i * self.masses[:, np.newaxis], axis=0) / self.masses.sum()
        eta_m = np.sum(eta_i * self.masses[:, np.newaxis], axis=0) / self.masses.sum()
        theta_m = np.arctan2(-eta_m, -xi_m) + np.pi
        # return self.wcoords.sum(axis=0) / self.masses.sum()

        return self.L * theta_m / 2 / np.pi

    @property
    def matrix(self):
        r"""
        Compute the gyration tensor. The diagonal terms are:

        1 / N \sum_i (x_i - x_com)^2

        for x, y and z, with com meaning center of mass and N the number of
        atoms.

        The off diagonal terms are:

        1 / N \sum_i (x_i - x_com) (y_i - y_com)
        """

        # mass weighted deviation with respect to the center of mass
        # x_i - x_com
        deviation = self.coords - self.center_of_mass
        # \sum_i (x_i - x_com)^2
        deviation2 = np.sum(deviation**2, axis=0)

        # off diagonal terms
        xy = np.dot(deviation[:, 0], deviation[:, 1])
        xz = np.dot(deviation[:, 0], deviation[:, 2])
        yz = np.dot(deviation[:, 1], deviation[:, 2])

        S = np.array([[deviation2[0], xy,            xz],
                      [xy,            deviation2[1], yz],
                      [xz,            yz,            deviation2[2]]])
        S /= self.masses.size

        return S

    @property
    def iso_rg(self):
        r"""
        Compute the radius of gyration, R_g, assuming an isotrop system as:

        R_g = \sqrt{1 / N \sum_i (r_i - r_com)^2}
        """
        natom = self.masses.size
        Rg2 = 1 / natom * np.sum((self.coords - self.center_of_mass)**2)

        return np.sqrt(Rg2)

    @property
    def iso_w_rg(self):
        r"""
        Compute the radius of gyration, R_g, assuming an isotrop system as:

        R_g = \sqrt{1 / N \sum_i (r_i - r_com)^2}
        """
        dr2 = (self.coords - self.center_of_mass)**2 * self.masses[:, np.newaxis]
        Rg2 = np.sum(dr2) / self.masses.sum()

        return np.sqrt(Rg2)

    @property
    def principal_moments(self):
        """
        Returns the sorted eigenvalue of the gyration tensor in decreasing
        order.
        """
        eigval, _ = np.linalg.eig(self.matrix)
        eigval = np.flip(np.sort(eigval))
        return eigval

    @property
    def rg(self):
        r"""
        Compute the radius of gyration from the eigenvalues of the gyration
        tensor from:

        R_g = \sqrt{\lambda_1 + \lambda_2 + \lambda_3}

        where the \lambda_i are the eigenvaleus of the gyration tensor.
        """
        return np.sqrt(np.sum(self.principal_moments))

    @property
    def shape_anisotropy(self):
        """
        Compute the shape anisotropy from:

        k2 = 1 - 3 (lambda_1 lambda_2 + lambda2 lambda3 + lambda3 lambda1) / (lambda1 + lambda2 + lambda3)^2
        """
        p1, p2, p3 = self.principal_moments
        return 1 - 3 * (p1 * p2 + p2 * p3 + p3 * p1) / (p1 + p2 + p3)**2

    @property
    def volume(self):
        r"""
        Compute the volume of an ellipsoid that would have its principal
        moments a, b and c equals to \sqrt(5 \lambda_i).
        """
        return 4 / 3 * np.pi * np.sqrt(5 ** 3 * self.principal_moments.prod())

    @property
    def density(self):
        """
        Compute the density assuming an effective ellipsoidal volume. The volume
        is assumed to be in nm^3 and the mass in g.mol-1. Thus the density is
        computed in kg.m-3
        """
        rho = np.sum(self.masses) / self.volume / avogadro / 1e-27 * 1e-3
        return rho

    @property
    def total_mass(self):
        """
        Compute the mass of aggregate in g.mol-1.
        """
        tmass = np.sum(self.masses)
        return tmass

    @property
    def max_distance(self):
        """
        Return the longest distance between an atom and the center of mass.
        """
        # Symmetrical distance matrix
        m = cdist(self.coords, self.coords, 'euclidean')

        # upper triangle matrix
        m = np.triu(m)
        # distances = np.sum((self.coords - self.center_of_mass)**2, axis=1)
        # return np.sqrt(np.max(distances))
        return np.max(m)

    @staticmethod
    def get_data_header():
        """
        Return an header according to the data returned by get_data
        """
        lines = "# Gyration tensor calculations.\n"
        lines += "# distance units depends on the units of the trajectory.\n"
        lines += "# (L) is the distance unit\n"
        lines += "# column 1: Rg ([L])\n"
        lines += "# column 2: Rg isotrop ([L])\n"
        lines += "# column 3: Rg isotrop mass weighted ([L])\n"
        lines += "# column 4: k^2, shape anisotropy (au) \n"
        lines += "# column 5: volume ([L]^3)\n"
        lines += "# column 5: molar mass ([g/mol])\n"
        lines += "# column 6: largest distance between center of mass and atoms ([L])\n"
        lines += "# column 7: density (g.L-1) grams per litters if [L] = nm.\n"
        # lines += "# column 8: Dipolar moment [D].\n"
        return lines

    def get_data(self):
        """
        Returns all data computed from the gyration tensor with 10.4f format.
        """
        lines = f"{self.rg:10.4f}"
        lines += f"{self.iso_rg:10.4f}"
        lines += f"{self.iso_w_rg:10.4f}"
        lines += f"{self.shape_anisotropy:10.4f}"
        lines += f"{self.volume:10.4f}"
        lines += f"{self.total_mass:12.4f}"
        lines += f"{self.max_distance:10.4f}"
        lines += f"{self.density:12.4e}"
        # lines += f"{self.mu:12.4e}"

        return lines


def center_of_mass(coords, L, masses=None):
    r"""Compute the center of mass of the points at coordinates `coords` with
    masses `masses`.
    Args:
        coords (np.ndarray): (N, 3) matrix of the points in :math:`\mathbb{R}^3`
        masses (np.ndarray): vector of length N with the masses
    Returns:
        The center of mass as a vector in :math:`\mathbb{R}^3`
    """
    # check coord array

    try:
        coords = np.array(coords, dtype=np.float64)
        coords = coords.reshape(coords.size // 3, 3)
    except ValueError:
        print("coords = ", coords)
        raise ValueError("Cannot convert coords in a numpy array of floats"
                         " with a shape (N, 3).")

    # check masses
    if masses is None:
        masses = np.ones(coords.shape[0])
    else:
        try:
            masses = np.array(masses, dtype=np.float64)
            masses = masses.reshape(coords.shape[0])
        except ValueError:
            print("masses = ", masses)
            raise ValueError("Cannot convert masses in a numpy array of "
                             "floats with length coords.shape[0].")
    if masses is None:
        masses = np.ones(coords.shape[0])

    """ compute the center of mass, the mass weighterd barycenter """

    theta_i = coords / L * 2 * np.pi
    xi_i = np.cos(theta_i)
    eta_i = np.sin(theta_i)
    xi_m = np.sum(xi_i * masses[:, np.newaxis], axis=0) / masses.sum()
    eta_m = np.sum(eta_i * masses[:, np.newaxis], axis=0) / masses.sum()
    theta_m = np.arctan2(-eta_m, -xi_m) + np.pi
    # return self.wcoords.sum(axis=0) / self.masses.sum()

    return L * theta_m / 2 / np.pi


def get_plane(coords, L, masses=None):
    r"""Given a set of N points in :math:`\mathbb{R}^3`, compute an orthonormal
    basis of vectors, the first two belonging to the plane and the third one
    being normal to the plane. In the particular case where N equal 3, there is
    an exact definition of the plane as the three points define an unique plan.
    If N = 3, use a gram-schmidt orthonormalization to compute the vectors. If
    N > 3, the orthonormal basis is obtained from SVD.
    Args:
        coords (np.ndarray): (N, 3) matrix of the points in :math:`\mathbb{R}^3`
        masses (np.ndarray): vector of length N with the masses
    Returns:
        Returns the orthonormal basis (vecx, vecy, n_a), vector n_a being
        normal to the plane.
    """
    # check coord array
    try:
        coords = np.array(coords, dtype=np.float64)
        coords = coords.reshape(coords.size // 3, 3)
    except ValueError:
        print("coords = ", coords)
        raise ValueError("Cannot convert coords in a numpy array of floats"
                         " with a shape (N, 3).")

    # check masses
    if masses is None:
        masses = np.ones(coords.shape[0])
    else:
        try:
            masses = np.array(masses, dtype=np.float64)
            masses = masses.reshape(coords.shape[0])
        except ValueError:
            print("masses = ", masses)
            raise ValueError("Cannot convert masses in a numpy array of "
                             "floats with length coords.shape[0].")

    com = center_of_mass(coords, L, masses)

    if coords.shape == (3, 3):
        # the plane is exactly defined from 3 points
        vecx = coords[1] - coords[0]
        vecx /= np.linalg.norm(vecx)

        # vecy, orthonormal with vecx
        vecy = coords[2] - coords[0]
        vecy -= np.dot(vecy, vecx) * vecx
        vecy /= np.linalg.norm(vecy)

        # normal vector
        n_a = np.cross(vecx, vecy)

    else:
        # get the best fitting plane from SVD.
        _, _, (vecx, vecy, n_a) = np.linalg.svd(coords - com)

    return vecx, vecy, n_a


def compute_dist_btw_com(coord_i, coord_j, L):

    theta_i = (coord_i / L) * 2 * np.pi
    theta_j = (coord_j / L) * 2 * np.pi

    xi_i = np.cos(theta_i)
    eta_i = np.sin(theta_i)

    xi_j = np.cos(theta_j)
    eta_j = np.sin(theta_j)

    xi_i_m = np.mean(xi_i, axis=0)
    eta_i_m = np.mean(eta_i, axis=0)

    xi_j_m = np.mean(xi_j, axis=0)
    eta_j_m = np.mean(eta_j, axis=0)

    theta_i_m = np.arctan2(-eta_i_m, -xi_i_m) + np.pi
    theta_j_m = np.arctan2(-eta_j_m, -xi_j_m) + np.pi

    # phi = np.abs(theta_i_m - theta_j_m)

    # compute norm
    phi = np.linalg.norm(theta_i_m - theta_j_m)

    # Use minimal image convention (in angles)

    # for degree more that 360
    if phi > 2 * np.pi:
        while phi > 2 * np.pi:
            phi -= 2 * np.pi

    if phi > np.pi:
        phi = (2 * np.pi) - phi

    if phi > np.pi / 2:
        phi = np.pi - phi

    return L * phi / 2 / np.pi


def angle_btw_flats(trajectory, imodel, jmodel, frame):
    """
    Return the distance between center of masse.

    """

    # 0. Selecting ariste box, assuming box cubic
    box = trajectory.unitcell_lengths[frame]
    L = box[0]

    # 1. selecting aromatic carbon.
    C_ar_i = trajectory.top.select(f"resid {imodel} and element C")
    C_ar_j = trajectory.top.select(f"resid {jmodel} and element C")

    # Selection coordinates
    coord_i = trajectory.xyz[frame][C_ar_i]
    coord_j = trajectory.xyz[frame][C_ar_j]

    # 2. getting the center of mass.
    com_i = center_of_mass(coord_i, L)
    com_j = center_of_mass(coord_j, L)

    # 3.1 calculing vector of plane
    _, _, n_i = get_plane(coord_i, L)
    _, _, n_j = get_plane(coord_j, L)

    # 3.2 calculing the distance between com-com
    #  d_com_com = np.linalg.norm(com_j - com_i)
    d_com_com = compute_dist_btw_com(coord_i, coord_j, L)

    # 4. angles between planes and angle between vector normal i and vector director from com_i to com_j
    u_n_i = n_i / np.linalg.norm(n_i)
    u_n_j = n_j / np.linalg.norm(n_j)

    u_com_ij = (com_j - com_i) / np.linalg.norm(com_j - com_i)

    angle_ni_nj = np.arccos(np.dot(u_n_i, u_n_j)) * 180.0 / np.pi

    ang_tp1 = np.arccos(np.dot(u_n_i, u_com_ij)) * 180.0 / np.pi
    ang_tp2 = np.arccos(np.dot(u_n_j, u_com_ij)) * 180.0 / np.pi

    if angle_ni_nj > 90.0:
        angle_ni_nj = 180.0 - angle_ni_nj

    # Selecting contact percentage

    # make pairs list
    # pairs = trajectory.top.select_pairs(i_ndx, j_ndx)

    return angle_ni_nj, d_com_com, ang_tp1, ang_tp2


def read_xtc(**kwargs):
    """ read the trajectrory using mdtraj module from the xtc file """

    # set needed variables:
    trajfile = kwargs["trajfile"]
    top = kwargs["top"]
    interval = kwargs["interval"]
    begin = kwargs["begin"]
    end = kwargs["end"]

    print("Read trajectory file: ", trajfile)
    print("Read topology from:   ", top)

    t1 = time.time()
    print("Start reading at", time.strftime("%H:%M:%S", time.localtime(t1)))
    print("  first frame: ", begin)
    print("  last frame:  ", end if end > 0 else "last")
    print("  step:        ", interval, "\n")

    trajectory = md.load(trajfile, top=top)
    t2 = time.time() - t1
    print("Done in %.0f s" % t2)

    print("  # atoms:           ", trajectory.n_atoms)
    print("  # frames total:    ", trajectory.n_frames)
    if trajectory.n_frames > 1:
        trajectory = trajectory[begin: end: interval]
    print("  # frames selected: ", trajectory.n_frames, "\n")

    return trajectory


def get_clusters(trajectory, ndx, rcut=0.5, ncut=1, interval=5, nmodels=5, **kwargs):
    """
    Look for clusters of residues defined in the given index according to the
    rcut and ncut criteria. Two residues belong to the same cluster if at
    least ncut distances between atoms of the two residues are lower than rcut.
    The clusters are identified for each frame of the trajectory.

    Args:
        trajectory (mdtraj.Trajectory): The trajectory
        ndx (list): 2D list. Provides the atom indexes of each residues. The
            first index run over residues, the second over atom index of the
            residue.
        rcut (float): cut-off distance that determines if residues belong to the
            same cluster. Unit are the same as the trajectory.
        ncut (int): number of distances that must be lower than rcut in order to
            gather residues in the same cluster.

    Returns:
        A list in which each element corresponding to one frame of the
        simulation. Each element corresponds to a list of cluster data. The
        whole strucutre looks like:
            [
                [
                    {"nmol": # of residue in the cluster,
                    "imol": list of residue ids in the cluster,
                    "ndx": [list of atom index of residues in the cluster]},
                    {"nmol": ...,
                    "imol": ...,
                    "ndx": [...]},
                    ...,
                ],
                [
                    {"nmol": ..., "imol": ..., "ndx": [...]},
                    {"nmol": ..., "imol": ..., "ndx": [...]},
                    ...,
                ],
                ...,
            ]
    """

    # compute distances between atoms in residue pairs for all frames
    print("Compute distances between residue pairs", end=" - ")
    t0 = time.time()

    n_models = len(ndx)
    distances = dict()

    for imodel, jmodel in it.combinations(range(n_models), 2):
        i_ndx = ndx[imodel]
        j_ndx = ndx[jmodel]

        # make pairs list
        pairs = trajectory.top.select_pairs(i_ndx, j_ndx)

        # compute distances with periodic boundary conditions
        ij_distances = md.compute_distances(trajectory, pairs)
        distances[(imodel, jmodel)] = ij_distances

    """
    for imodel in range(n_models - 1):
        i_ndx = ndx[imodel]
        for jmodel in range(imodel + 1, n_models):
            j_ndx = ndx[jmodel]

            # make pairs list
            pairs = trajectory.top.select_pairs(i_ndx, j_ndx)

            # compute distances with periodic boundary conditions
            ij_distances = md.compute_distances(trajectory, pairs)
            distances[(imodel, jmodel)] = ij_distances
    """

    dt = time.time() - t0
    print("Done in %.0f s\n" % dt)

    # look for the number of clusters
    print("Start cluster analysis")
    print(f"  rcut: {rcut:8.3f}")
    print(f"  ncut: {ncut:8d}")
    print(f"  nmol: {nmodels:8d}")

    # list to store cluster of each frame
    clusters = list()

    # dictionary with angles history between two planes in molecules
    # the plane are taked from carbon aromatic from core.
    # angles_pairs = [iframe distance mol1 mol2 ang1 dist_com ang_mol1_com ang_mol2_com ]
    angles_pairs = list()

    t1 = time.time()
    t0 = t1
    for iframe in range(trajectory.n_frames):
        # time evaluation
        true_frame = iframe * interval
        if iframe != 0 and iframe % (trajectory.n_frames // 10) == 0:
            t2 = time.time()
            dt = t2 - t1
            percentage = iframe / trajectory.n_frames * 100
            eta = dt * 10 - (t2 - t0)
            print(f" **** {percentage:3.0f}%% in {dt:.0f} s - "
                  f"frame {iframe:4d} (ETA {eta:.0f} s)")
            t1 = t2

        frame_clusters = list()
        for (imodel, jmodel) in distances:
            pair_distances = distances[(imodel, jmodel)]

            # look for how many distances are lower than rcut value
            n_neighbors = np.where(pair_distances[iframe] < rcut)[0].size

            if n_neighbors >= ncut:
                # at least ncut distances are lower than rcut
                # imodel and jmodel belong to the same cluster
                # The number of neighbors needed to be in the same cluster can be
                # modified

                n_at_tot = np.sqrt(len(pair_distances[iframe]))
                p_contact = n_neighbors * 100 / n_at_tot

                # initialize clusters
                if frame_clusters == []:
                    frame_clusters.append([imodel, jmodel])

                    d_ij_min = distances[(imodel, jmodel)][iframe].min()
                    ang_ni_nj, d_com_com, ang_tp1, ang_tp2 = angle_btw_flats(trajectory, imodel, jmodel, iframe)

                    # [iframe, distance_min, mol1, mol2, ang_ni_nj, dist_com, ang_mol1_com, ang_mol2_com ]
                    angles_pairs.append([true_frame, d_ij_min, imodel, jmodel, ang_ni_nj, d_com_com, p_contact])

                else:
                    # look for the cluster containing imodel and jmodel
                    icluster = None
                    jcluster = None
                    for ic, cluster in enumerate(frame_clusters):
                        if imodel in cluster:
                            icluster = ic
                        if jmodel in cluster:
                            jcluster = ic

                    # print("\nModel: ", (imodel, jmodel), (icluster, jcluster), end= " - ")
                    # update clusters list
                    if (icluster is None) and (jcluster is None):
                        # create a new cluster
                        # print("new cluster")
                        frame_clusters.append([imodel, jmodel])

                        d_ij_min = distances[(imodel, jmodel)][iframe].min()
                        ang_ni_nj, d_com_com, ang_tp1, ang_tp2 = angle_btw_flats(trajectory, imodel, jmodel, iframe)

                        # [iframe, distance_min, mol1, mol2, ang_ni_nj, dist_com, ang_mol1_com, ang_mol2_com ]
                        angles_pairs.append([
                            int(true_frame), float(d_ij_min), int(imodel), int(jmodel),
                            float(ang_ni_nj), float(d_com_com), p_contact])

                    elif (icluster is not None) and (jcluster is None):
                        # add jmodel to the cluster of imodel
                        # print(f"add jmodel {jmodel} in the cluster {icluster} with imodel {imodel}")
                        frame_clusters[icluster].append(jmodel)

                        d_ij_min = distances[(imodel, jmodel)][iframe].min()
                        ang_ni_nj, d_com_com, ang_tp1, ang_tp2 = angle_btw_flats(trajectory, imodel, jmodel, iframe)

                        # [iframe, distance_min, mol1, mol2, ang_ni_nj, dist_com, ang_mol1_com, ang_mol2_com ]
                        angles_pairs.append([
                            int(true_frame), float(d_ij_min), int(imodel), int(jmodel),
                            float(ang_ni_nj), float(d_com_com), p_contact])

                    elif (icluster is None) and (jcluster is not None):
                        # add imodel to the cluster of jmodel
                        # print(f"add imodel {imodel} in the cluster {jcluster} with jmodel {jmodel}")
                        frame_clusters[jcluster].append(imodel)

                        d_ij_min = distances[(imodel, jmodel)][iframe].min()
                        ang_ni_nj, d_com_com, ang_tp1, ang_tp2 = angle_btw_flats(trajectory, imodel, jmodel, iframe)

                        # [iframe, distance_min, mol1, mol2, ang_ni_nj, dist_com, ang_mol1_com, ang_mol2_com ]
                        angles_pairs.append([
                            int(true_frame), float(d_ij_min), int(imodel), int(jmodel),
                            float(ang_ni_nj), float(d_com_com), p_contact])

                    # at this step neither icluster nor jcluster equal None
                    elif icluster != jcluster:
                        # merge clusters
                        # print(f"merge clusters {icluster} and {jcluster}")
                        merged = frame_clusters[icluster] + frame_clusters[jcluster]
                        frame_clusters = [fc for ic, fc in enumerate(frame_clusters)
                                          if ic not in (icluster, jcluster)]
                        frame_clusters += [merged]

        # complete frame_clusters with isolated models
        models_in_clusters = set([ic for cluster in frame_clusters for ic in cluster])
        for ic in range(n_models):
            if ic not in models_in_clusters:
                frame_clusters.append([ic])

        # set up clusters data
        data = list()
        for cluster in frame_clusters:
            data.append(
                dict(nmol=len(cluster), imol=cluster,
                     ndx=[ndx[imodel] for imodel in cluster])
            )

        clusters.append(data)

    dt = time.time() - t0
    print("Done in %.0f s\n" % dt)

    return clusters, angles_pairs


def get_options():
    """ get options from the command line """

    parser = argparse.ArgumentParser(
        prog="gro2gau",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    # trying to mimic gromacs gmx option names
    parser.add_argument("-f", "--trajfile",
                        help="GROMACS trajectory file.",
                        metavar="TRAJFILE",
                        default="traj_comp.xtc",
                        type=str)

    parser.add_argument("-s", "--top",
                        help="A file which contains topology: a gro or pdb file for example.",
                        metavar="TOP",
                        default="confout.gro",
                        type=str)

    parser.add_argument("-i", "--interval",
                        help="Interval between frames. (default = 1)",
                        metavar="I",
                        default=1,
                        type=int)

    parser.add_argument("-b", "--begin",
                        help="Index of the first frame. (default 0)",
                        metavar="B",
                        default=0,
                        type=int)

    parser.add_argument("-e", "--end",
                        help="Index of the last frame. (default last)",
                        metavar="E",
                        default=-1,
                        type=int)

    parser.add_argument("-r", "--rcut",
                        help=("Cut-off distance. Two residues belong to the same "
                              "cluster if at least ncut distances are lower than "
                              "rcut. Same units as the trajectory."),
                        metavar="RCUT",
                        default=0.5,
                        type=float)

    parser.add_argument("-n", "--ncut",
                        help=("Number of distances between atoms of two "
                              "residues that must be lower than rcut to decide "
                              "if these two residues belong to the same cluster."),
                        metavar="NCUT",
                        default=1,
                        type=int)

    parser.add_argument("-nmol", "--nmodels",
                        help=("Number of models selected."),
                        metavar="NMOL",
                        default=20,
                        type=int)

    return vars(parser.parse_args())


if __name__ == "__main__":

    # starting setup
    args = get_options()
    print(TITLE)

    # read trajectory
    trajectory = read_xtc(**args)

    # make index for models
    # this is a crucial point to reduce computational cost
    n_models = args["nmodels"]
    models_ndx = np.array(
        [trajectory.top.select(f"resid {ires} and not name HC") for ires in range(n_models)]
    )

    # compute clusters data
    clusters, angles_pairs = get_clusters(trajectory, models_ndx, **args)
    masses = np.array([atom.element.mass for atom in trajectory.topology.atoms])

    # file with gyration data
    lines = GyrationTensor.get_data_header()
    lines += "# column 8: number of models in the cluster\n"
    lines += "# column 9: frame number\n"

    # file with the number of cluster per frame
    nclust_lines = "# Number of clusters in each frame\n"

    # file with the residue id for each cluster and each frame
    resid_lines = "# Residue id for each frame and for each cluster\n"
    resid_lines += "# iframe   nmol   [resid]\n"

    print("Analysis evolution:")
    t1 = time.time()
    t0 = t1
    for iframe, frame_clusters in enumerate(clusters):
        true_frame = iframe * args["interval"]

        if iframe != 0 and iframe % (trajectory.n_frames // 10) == 0:
            t2 = time.time()
            dt = t2 - t1
            percentage = iframe / trajectory.n_frames * 100
            eta = dt * 10 - (t2 - t0)
            print(f" **** {percentage:3.0f}%% in {dt:.0f} s - "
                  f"frame {iframe:4d} (ETA {eta:.0f} s)")
            t1 = t2

        try:
            nclust_lines += f"{true_frame:8d} {len(frame_clusters):5d}\n"

            for cluster in frame_clusters:
                ndx = np.hstack(cluster["ndx"])

                gyr = GyrationTensor(
                    coords=trajectory.xyz[iframe, ndx],
                    masses=masses[ndx],
                    box=trajectory.unitcell_lengths)
                lines += gyr.get_data()

                lines += f"{cluster['nmol']:4d}"
                lines += f"{true_frame:8d}\n"

                resid_lines += f"{true_frame:8d}{cluster['nmol']:4d}"
                resid_lines += "".join([f"{resid:4d}" for resid in cluster["imol"]])
                resid_lines += "\n"

        except TypeError:
            print(true_frame)
            continue

    angles_lines = "# Data for angles between pairs of molecules.\n"
    angles_lines += "# iframe   d_ij_min   imol   jmol   angle_ni_nj   d_com_com   p_contact\n"

    for pairs in angles_pairs:
        angles_lines += f"{pairs[0]:8d}"
        angles_lines += f"{pairs[1]:11.4f}"
        angles_lines += f"{pairs[2]:7d}"
        angles_lines += f"{pairs[3]:7d}"
        angles_lines += f"{pairs[4]:14.2f}"
        angles_lines += f"{pairs[5]:12.4f}"
        angles_lines += f"{pairs[6]:12.2f}"
        angles_lines += "\n"

    print("Save data into files.")

    # write files
    with open("clusters.dat", "w") as f:
        f.write(lines)

    with open("nclusters.dat", "w") as f:
        f.write(nclust_lines)

    with open("cluster_resid.dat", "w") as f:
        f.write(resid_lines)

    with open("angles_pairs.dat", "w") as f:
        f.write(angles_lines)
