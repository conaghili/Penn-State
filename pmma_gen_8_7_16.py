#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import permutations as perm

"""
pmma_gen.py

This script generates a coarse-grained model of poly-methylmethacrylate (PMMA)
compatible with the NAMD molecular dynamics software package.

The user can set how many coarse-grained PMMA chains are generated, how many
monomers per chain, as well as the simulation cell size.

Files that are generated:
    - cg_pmma.pdb
    - cg_pmma.top
"""

def main():
    """ main function used to generate inputs """

    # -- Simulation cell parameters
    number_of_chains = 20
    chain_length = 75  # -- monomers
    box_length = 10    # -- in units of [nm]

    # -- Model generation -- Do not modify below unless you are sure --
    global bead_id, chain_id
    bead_id, chain_id = 1, 1

    polymer_system = build_chains(chain_length, number_of_chains, box_length)
    write_pdb(chain_length, polymer_system)
    write_fixed_pdb(chain_length, polymer_system)
    write_top_file(chain_length, polymer_system)
    

    return

def write_pdb(chain_length, polymer_system):
    """
    Save bead positions to a pdb format

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)
    [2] http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    """

    with open("index.pdb", "w") as out_file:
        out_file.write(SEQRES(polymer_system, chain_length))
                 
        for chain in polymer_system:
            chainIDval = chain[0].chain_id
            resSeq = 1
            for bead in chain:
               
                # -- Define parameters to describe a bead in pdb file
                chain_label = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                               'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                               'U', 'V', 'W', 'X', 'Y', 'Z']

                chainID = chain_label[chainIDval-1]
                x = bead.position[0]
                y = bead.position[1]
                z = bead.position[2]
                Occupancy = 1.0
                TempFactor = 0.0
                bead.id_num = bead.id_num  

                # -- The chain ID is typically a one string character, like 'X'
                # -- We're trying to use an int for the chain id for the moment
                count = 3000
                out_file.write('{:<6}{:>5d}{:>3}{:>6}{:>2}{:>4}{:12.3f}' \
                    '{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}{:>2}\n'.format("ATOM",
                    bead.id_num, bead.cg_type, "MON", chainID, resSeq, x, y,
                    z, Occupancy, TempFactor, bead.cg_type, "0"))
                if bead.cg_type == 'B': resSeq += 1

       
        
          # -- Define the structural parameters of an FCC lattice
        a = 0.15*np.sqrt(2.)
        a1 = 0.5*a*np.array([1., 1., 0.])
        a2 = 0.5*a*np.array([1., 0., 1.])
        a3 = 0.5*a*np.array([0., 1., 1.])
        basis = np.array([a1, a2, a3])

        # -- Generate a lattice
        x = 3
        perm_range = range(-x,x) + range(-x,x) + range(-x,x)
        points = set( perm( perm_range, 3 ) )
        points = np.array( list(points) )
        lattice = np.matmul(points, basis)

        # -- Cut out a sphere from the lattice
        radius = 0.3
        new_lattice = np.array([])
        for point in lattice:
            if np.linalg.norm(point) < radius:
                new_lattice = np.append(new_lattice, point)

        new_lattice = new_lattice.reshape(new_lattice.size/3, 3)
        
        # -- Plot the sphere
        #norms = np.linalg.norm(new_lattice, axis=1)
    
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')

        #ax.scatter(new_lattice[:,0], new_lattice[:,1], new_lattice[:,2], c=norms, s=200)
        #ax.set_aspect('equal')
        #plt.show()

        for row in new_lattice: # -- each row in the lattice stores the coordinate for an atom
            count = count  +1
            out_file.write('{:<6}{:>5d}{:>3}{:>6}{:>2}{:>4}{:12.3f}' \
                    '{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}{:>2}\n'.format("ATOM",
                    count, "C", "NAN","U",count, row[0] +4.5, row[1] +4.5, row[2] +4.5,Occupancy, 
                    TempFactor,"C", "0"))
            
        out_file.write( '{:<6}'.format('END') )
        
""" 
fixed pdb file
 """
                
def write_fixed_pdb(chain_length, polymer_system):        
 with open("index_fixed.pdb", "w") as out_file:
        out_file.write("""
""")
        for chain in polymer_system:
            chainIDval = chain[0].chain_id
            resSeq = 1
            for bead in chain:
               
                # -- Define parameters to describe a bead in pdb file
                chain_label = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                               'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                               'U', 'V', 'W', 'X', 'Y', 'Z']

                chainID = chain_label[chainIDval-1]
                x = bead.position[0]
                y = bead.position[1]
                z = bead.position[2]
                Occupancy = 1.0
                TempFactor = 0.0
                bead.id_num = bead.id_num 
                count = 3000

                # -- The chain ID is typically a one string character, like 'X'
                # -- We're trying to use an int for the chain id for the moment

                out_file.write('{:<6}{:>5d}{:>3}{:>6}{:>2}{:>4}{:12.3f}' \
                    '{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}{:>2}\n'.format("ATOM",
                    bead.id_num, bead.cg_type, "MON", chainID, resSeq, x, y,
                    z, Occupancy, TempFactor, bead.cg_type, "0"))
                if bead.cg_type == 'B': resSeq += 1

          
          # -- Define the structural parameters of an FCC lattice
        a = 1*np.sqrt(2.)
        a1 = 0.5*a*np.array([1., 1., 0.])
        a2 = 0.5*a*np.array([1., 0., 1.])
        a3 = 0.5*a*np.array([0., 1., 1.])
        basis = np.array([a1, a2, a3])

        # -- Generate a lattice
        x = 8
        perm_range = range(-x,x) + range(-x,x) + range(-x,x)
        points = set( perm( perm_range, 3 ) )
        points = np.array( list(points) )
        lattice = np.matmul(points, basis)

        # -- Cut out a sphere from the lattice
        radius = 6.
        new_lattice = np.array([])
        for point in lattice:
            if np.linalg.norm(point) < radius:
                new_lattice = np.append(new_lattice, point)

        new_lattice = new_lattice.reshape(new_lattice.size/3, 3)
        
        # -- Plot the sphere
        #norms = np.linalg.norm(new_lattice, axis=1)
    
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')

        #ax.scatter(new_lattice[:,0], new_lattice[:,1], new_lattice[:,2], c=norms, s=200)
        #ax.set_aspect('equal')
        #plt.show()

        for row in new_lattice: # -- each row in the lattice stores the coordinate for an atom
            count = count  +1
            Temp = 1
            out_file.write('{:<6}{:>5d}{:>3}{:>6}{:>2}{:>4}{:12.3f}' \
                    '{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}{:>2}\n'.format("ATOM",
                    count, "C", "NAN","U",count, row[0] + 4.5, row[1] + 4.5, row[2] + 4.5,Occupancy, 
                    Temp,"C", "0"))
            
        
        out_file.write( '{:<6}'.format('END') )
        
        


def write_top_file(chain_length, polymer_system):
    """
    Write the topology file

    References
    ----------
    [1] https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/rtop/
    """


    topology = """
    *>>>>>> CHARMM Topology File for Coarse-Grained PMMA <<<<<<<<<
*
27  1

! References:
!
! CHARMM Topology file format
! ----------------------------
! https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/rtop/
!
!

MASS     1 C1      72.0000 A ! carbon backbone bead
MASS     2 Na      72.0000 B ! R-group bead
MASS     3 P1      72.0000 C ! Nano Particle 
!  DECL declares which atoms are selected from previous or next residues to
!  join residues together
DECL -A
DECL +A
DEFA FIRS NTER LAST CTER

RESI MON         0.00
GROUP
ATOM A    C1     0.00  !     ( A – B )_n
ATOM B    Na     0.00  !
BOND A B  A +A

RESI NAN         0.00
GROUP
ATOM C    P1     0.00  !    
END
    """
    with open("pmma_CG_topology.inp", 'w') as out:
        out.write(topology)


def SEQRES(polymer_system, chain_length):
    """
    returns a formatted string containing the residue sequence block

    13 residues per line

             1         2         3         4         5         6         7         8
    12345678901234567890123456789012345678901234567890123456789012345678901234567890
    SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU
    SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN
    SEQRES   1 B   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU
    SEQRES   2 B   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR
    SEQRES   3 B   30  THR PRO LYS ALA
    SEQRES   1 C   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU
    SEQRES   2 C   21  TYR GLN LEU GLU ASN TYR CYS ASN
    SEQRES   1 D   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU
    SEQRES   2 D   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR                                                                                                                                SEQRES   3 D   30   THR PRO LYS ALA
    SEQRES   1 A    8   DA  DA  DC  DC  DG  DG  DT  DT
    SEQRES   1 B    8   DA  DA  DC  DC  DG  DG  DT  DT
    SEQRES   1 X   39    U   C   C   C   C   C   G   U   G   C   C   C   A
    SEQRES   2 X   39    U   A   G   C   G   G   C   G   U   G   G   A   A
    SEQRES   3 X   39    C   C   A   C   C   C   G   U   U   C   C   C   A
    """

    serNum = 1
    chainID = 1
    # -- Not sure if we are limited to single characters for the chain length
    chain_label = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J','K', 'L',
                   'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
                   'X', 'Y', 'Z']
    sequence = ''
    for chain in polymer_system:
        chainID = chain[0].chain_id
        number_of_records = int(chain_length // 13. + 1)
        remaining_monomors = chain_length - (number_of_records - 1)*13
        for serNum in range(1,number_of_records+1):
            if number_of_records != 1 and serNum < number_of_records:
                record="{:<9}{:<2}{:<4}{:<4}".format("SEQRES",serNum,
                                             chain_label[chainID-1],chain_length)
                temp = str()
                seq = [ temp + "{:<4}".format("MON") for i in range(13) ]
                seq = ''.join(seq)+'\n'
                record += seq
            if serNum == number_of_records:
                record="{:<9}{:<2}{:<4}{:<4}".format("SEQRES",serNum,
                                             chain_label[chainID-1],chain_length)
                temp = str()
                seq = [temp + "{:<4}".format("MON") for i in range(remaining_monomors)]
                seq = ''.join(seq)+'\n'
                record += seq
            sequence += record
    return sequence




def build_chains(chain_length, number_of_chains, box_length):
    """
    Parameters
    ----------
    chain_length : int
        Number of A-B pairs in the chain

    box_length : float
        Length of each side of the simulation box

    Returns
    -------
    polymer_system :
        List of beads

    """

    global bead_id, chain_id

    # -- Initialize a list to contain the polymer system
    system = list()

    # -- Now build the chains
    for chain_id in range(1, number_of_chains+1):

        chain = list()

        # -- Build chain one monomer at a time
        for i in range(1, chain_length+1):
            print "i = %d" % (i)
            if i is 1:
                # -- Calculate the first 'A' position of the chain randomly
                a_position = np.random.rand(3)*box_length
                chain.append(CGBead(bead_id, chain_id, 'A', a_position))

                # -- The first 'B' bead is also a random vector
                bead_id += 1
                b_position = a_position + random_vect()
                chain.append(CGBead(bead_id, chain_id, 'B', b_position ))

            elif i is 2:
                # -- The next 'D' bead placement is constrained by the D-D-B angle
                temp = position_n1b - position_n1a
                a_position = position_n1a + get_next_position(temp, 71.)
                chain.append(CGBead(bead_id, chain_id, 'A', a_position))

                # -- The next 'B' bead placement is constrained by the D-D-B angle
                bead_id += 1
                temp = position_n1a - a_position
                b_position = a_position + get_next_position(temp, 71.)
                chain.append(CGBead(bead_id, chain_id, 'B', b_position))

            else:
                # -- The next 'A' bead placement is constrained by both the D-D-D
                # -- angle and the A-A-B angle
                temp_1 = position_n2a - position_n1a
                temp_2 = position_n1b - position_n1a
                #theta = 38.41786204357107
                a_position = get_next_position(temp_1, 131.) + position_n1a
                chain.append(CGBead(bead_id, chain_id, 'A', a_position))

                # -- The next 'B' bead placement is constrained by the D-D-B angle
                bead_id += 1
                temp = position_n1a - a_position
                b_position = a_position + get_next_position(temp, 71.)
                chain.append(CGBead(bead_id, chain_id, 'B', b_position))

            bead_id += 1
            position_n1a = chain[2*i-2].position
            position_n1b = chain[2*i-1].position
            print "a_position = " + str( position_n1a )
            print "b_position = " + str( position_n1b )
            print ""
            if i > 1:
                position_n2a = chain[2*i-4].position
                position_n2b = chain[2*i-3].position

        system += [chain]
    return system

def get_next_position(a, theta, constrained = False, b = None):
    """
    This function calculates a random vector that satisfies the
    bond angle requirements for the backbone A beads (131 degrees).

    Parameters
    ----------
    a : np.array, shape (3,)
        Vector describing position between previous two beads

    theta : float
        Angle in degrees

    Returns
    -------
    x : np.array, shape (3,)
        Vector for next bead position
    """

    # -- Define the basic structural parameters of the beads
    r = 0.282
    theta *= np.pi/180.  # -- convert to radians
    phi = np.random.rand()*2*np.pi
    if constrained:
        theta_ac = 131. * np.pi / 180.
        theta_bc = 71. * np.pi / 180.
        phi = compute_phi(r, b, theta_ac, theta_bc)

    # -- Compute the point rotated w.r.t. the z-axis (from spherical coord.)
    X = r*np.array([np.sin(theta)*np.cos(phi),
                    np.sin(theta)*np.sin(phi),
                    np.cos(theta)])

    # -- Rotate the bead coord. so that the z-axis is aligned with 'a'
    # -- This ensures the angle is always correct and position random
    Z = np.array([0., 0., 1.])
    a /= np.linalg.norm(a)
    rotation_axis = np.cross(a, Z)
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_angle = get_angle(a, Z)*-1
    X = rotate(X, rotation_axis, rotation_angle)
    # print get_angle(a,X)*180./np.pi
    return X

def compute_phi(r, b, theta_ac, theta_bc):
    """
    Reference : https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
    """
    arg = np.linalg.norm(b)**2 * np.cos(theta_bc) - b[2]*r*np.cos(theta_ac)
    arg *= 1. / r / np.sqrt(b[0]**2 + b[1]**2) / np.sin(theta_ac)
    print "arg = %f" % (arg)
    return np.arcsin(arg) - np.arctan2(b[0], b[1])

def random_vect():
    """
    Define a random vector of length 0.282 nm.
    This is the equilibrium length of any D-D or D-B bond.
    """
    rand_vector = np.random.rand(3)
    rand_vector /= np.linalg.norm(rand_vector)
    return rand_vector*0.282

def get_angle(a,b):
    """
    returns the angle in radians between vectors a and b
    """
    arg = np.dot(a,b) / np.linalg.norm(a) / np.linalg.norm(b)
    # -- Because of floating point precision, angles of n*pi are treated poorly
    if np.isclose(arg, 1): arg = 1
    if np.isclose(arg, -1): arg = -1
    return np.arccos(arg)

def rotate(points, rotation_axis, theta, angle_unit='radians'):
    """
    Rotates a set of vectors about an arbitrary axis of rotation
    and by an arbitrary angle.

    Parameters
    ----------

    points : np.ndarray, [N x 3]
        A set of points to be rotated

    rotation_axis : np.ndarray, [1 x 3]
        The axis <u,v,w> about which the points are rotated.

    theta : float
        The angle by which the points are rotated. Can be in
        units of degrees or radians.

    Kwargs
    ------

    angle_unit : string, optional
        Controls whether the angle is in units of degrees or radians.
        Options are : 'degrees', 'radians'

    Returns
    -------

    rotated_points : np.ndarray, [N x 3]
        The set of rotated points

    Reference(s)
    ------------

    1) http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    Notes
    -----

    We assume for the moment that the axis of rotation passes through
    the origin. The true general function would specify two points P1 and P2
    that the axis of rotation passes through so that the axis of
    rotation is given by: <u,v,w> = <P2x - P1x, P2y - P1y, P2z - P1z> .
    """

    X = points.T

    u = rotation_axis[0]
    v = rotation_axis[1]
    w = rotation_axis[2]

    # -- For the moment we assume the axis of rotation passes through the origin
    a = 0.0
    b = 0.0
    c = 0.0

    X_prime = np.zeros_like(X)
    X_prime[0] = (a*(v**2 + w**2) - u*(b*v + c*w - u*X[0] - v*X[1] - w*X[2]))* \
                 (1 - np.cos(theta)) + X[0]*np.cos(theta) + \
                 (-c*v + b*w - w*X[1] + v*X[2])*np.sin(theta)
    X_prime[1] = (b*(u**2 + w**2) - v*(a*u + c*w - u*X[0] - v*X[1] - w*X[2]))* \
                 (1 - np.cos(theta)) + X[1]*np.cos(theta) + \
                 (c*u - a*w + w*X[0] - u*X[2])*np.sin(theta)
    X_prime[2] = (c*(u**2 + v**2) - w*(a*u + b*v - u*X[0] - v*X[1] - w*X[2]))* \
                 (1 - np.cos(theta)) + X[2]*np.cos(theta) + \
                 (-b*u + a*v - v*X[0] + u*X[1])*np.sin(theta)

    return X_prime.T

class CGBead:
    """ Defines a bead object. """
    def __init__(self, bead_id, chain_id, bead_type, position):
        self.id_num = bead_id
        self.chain_id = chain_id
        self.cg_type = bead_type
        self.position = np.array(position)

if __name__ == "__main__":
    main()
