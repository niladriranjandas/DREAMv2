"""Helper module of <m torsionInterpolPot>.

Tools for setup of torsion angle interpolated potential terms.


When using this module, please cite:

Bermejo, G.A., Clore, G.M., and Schwieters, C.D. (2012). Smooth statistical
torsion angle potential derived from a large conformational database via
adaptive kernel density estimation improves the quality of NMR protein
structures. Protein Sci. 21, 1824-1836.

and possibly also:

Hu, K.N., Qiang, W., Bermejo, G.A., Schwieters, C.D., and Tycko, R. (2012).
Restraints on backbone conformations in solid state NMR studies of uniformly
labeled proteins from quantitative amide N-15-N-15 and carbonyl C-13-C-13
dipolar recoupling data. J. Magn. Reson. 218, 115-127.

"""
# Written by Guillermo A. Bermejo.

import os
import glob

import torsionInterpolPot  
import potList
import cdsVector
import simulation



def getdatafilepath(data=None, ext='.dat'):
    """Return a list with the path(s) of data file(s).

    Argument data is a string with the name of either a data file or a directory
    containing one or more data files (with filename extention given by the ext
    argument, a string). Alternatively, a sequence with several such filenames
    can be provided as the data argument.  If data is not provided, all files in
    the current directory with filename extention matching that ext are assumed
    to be data files.

    """ 
    if data is None:
        result = glob.glob('./*'+ext)
    else:
        if type(data) is str: data = [data]
        result = []
        for name in data:
            if os.path.isfile(name): result.append(name)
            if os.path.isdir(name): result.extend(glob.glob(name+'/*'+ext))
    if result:
        return result
    else:
        raise Exception('No data files found')



def setup_pot(name, info, data, selstrings=None, verbose=False, sim=None):
    """Return an N-dimensional <m torsionInterpolPot>.TorsionInterpolPotND.

    The name argument (a string) is the name of the potential.  info is a
    comment string that can be retrieved from the returned potential via the
    comment method.  The data argument is a dictionary that contains all the
    numerical data of the potential; its keys are integers used to identify the
    different torsion angles and energy values.  For example, in a 2D case we
    may have keys 1, 2, and 999, respectively corresponding to the first and
    second torsion angle, and the energy values (the largest key labels the
    energy).  The numbers associated with the torsion angle axes and energy are
    provided by the dict values in the form of a list.  The selstrings argument
    is a list of strings, each an XPLOR selection string for an atom that, along
    with the rest, define the torsion angle(s) to which the potential applies.
    For example,

        selstrings = ["name C and bondedto (resid 25 and name N)",
                      "resid 25 and name N",
                      "resid 25 and name CA",
                      "resid 25 and name C",     # last atom of first torsion
                      "resid 25 and name N",     # first atom of second torsion
                      "resid 25 and name CA",
                      "resid 25 and name C",
                      "name N and bondedto (resid 25 and name C)"]

    applies the 2D potential to phi and psi of residue 25.  The order of the
    torsion angles in selstrings corresponds to that given by the keys of the
    data dict.  If selstrings is not provided, the potential term is still
    returned, but it does not apply to any torsion.  If verbose is True, some
    information is printed about the setup process.  sim is a simulation object
    which defaults to <m simulation>.currentSimulation().
    
    """
    # Simulation.
    if sim is None: sim = simulation.currentSimulation()  
    
    # Generate list of numbers ("nums") from data dict.
    # For example, in a 2D case nums = [x, y, energy], where each element 
    # (e.g., x) is a list with the correnponding numbers.
    nums = list(data.items())
    nums.sort()
    nums = [t[1] for t in nums]
    
    if selstrings is not None:
        if len(selstrings)//4 != len(nums)-1:
            raise Exception('Mismatch between dimensionality from \
selstrings (%i) and torsion angle axes (%i)' % (len(selstrings)/4, len(nums)-1))

    ndim = len(nums) - 1  # number of dimensions

    # Create potential term.
    errmsg = 'Mismatch between the provided number of values in axes and energy'

    if ndim == 1:

        if len(nums[1]) != len(nums[0])-1:
            raise Exception(errmsg)
        # Directly inputting list versions to
        # torsionInterpolPot.TorsionInterpolPot1D fails.
        x = cdsVector.CDSVector_double(nums[0])
        energy = cdsVector.CDSVector_double(nums[1])
        pot = torsionInterpolPot.TorsionInterpolPot1D(name, x, energy, sim)
        if info: pot.setComment(info)
        
    elif ndim == 2:

        if len(nums[2]) != (len(nums[0])-1)*(len(nums[1])-1):
            raise Exception(errmsg)
        pot = torsionInterpolPot.TorsionInterpolPot2D(name, nums[0], nums[1],
                                                      nums[2], sim)
        if info: pot.setComment(info)

    elif ndim == 3:

        if len(nums[3]) != ((len(nums[0])-1)*(len(nums[1])-1)*(len(nums[2])-1)):
            raise Exception(errmsg)
        pot = torsionInterpolPot.TorsionInterpolPot3D(name, nums[0], nums[1],
                                                      nums[2], nums[3], sim)
        if info: pot.setComment(info)
        
    else:
        raise Exception('%i is an unsupported dimensionality' % ndim)

    # Potential term created. If non-empty selstrings is provided, add the
    # restraint to it (assumes the same interface for all dimensionalities).

    # Note that this function allows the addition of only one restraint to 
    # the potential.

    if selstrings:
        # For lack of better option, restraint has same name as the potential.
        restraint = 'name %s\nassign %s' % (name, ' '.join(selstrings))
        pot.addRestraints(restraint)  # will fail if any selstring selects
                                      # other than one atom
        if verbose: print('%sD torsionInterpolPot named %s created' % (ndim,
                                                                       name))
    return pot



def create_TorsionInterpolPot(name='', data=None, ext='.dat', verbose=False):
    """Return a <m potList>.PotList of torsion angle potential terms.

    Sets up torsion angle interpolated potential terms of the type defined in
    the <m torsionInterpolPot> module.  The name argument gives the name to the
    returned <m potList>.PotList.  data argument is a string with the name of
    either a data file containing the energy points to be interpolated, or a
    directory containing one or more such data files (in the latter case, only
    filenames with the extension given by the ext argument are considered).
    Alternatively, a sequence with several such filenames can be provided as the
    data argument.  For documentation on the format of data files see the
    README.format file in the databases/torsions directory within the Xplor-NIH
    directory.  If verbose=True, some information is printed about the setup
    process.

    Each data file (indicated by the data argument) may specify one or more
    potential terms, which are packed into a <m potList>.PotList. The latter is
    appended to the returned potlist, which thus consists of a potlist of
    potlists.

    """
    result = potList.PotList(name)
    filenames = getdatafilepath(data, ext)

    potnames = {} # to count the potentials

    for filename in filenames:

        if verbose: print('Processing', filename)

        # Create a potlist for all terms within this file.
        name = filename.split('/')[-1] # get filename without path
        name = name[:name.find('.')]  # discard filename extension
        pots = potList.PotList(name)  # to hold potential terms in current file
        pots.resetPotName('TorPotl')
        
        infile = open(filename, 'r')
        contents = infile.read()
        infile.close()

        name = ''  # now the name of the potential term (not the result potlist)
        lines = contents.splitlines()
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#'):
                if line.lstrip().startswith('name'):
                    if name: # the start of a new term; set up the previous one
                        pot = setup_pot(name, '\n'.join(info), data, selstrings,
                                        verbose)
                        pots.append(pot)
                        potname = pot.potName()
                        potnames[potname] = potnames.setdefault(potname, 0) + 1
                    name = line.split()[1]
                    selstrings = []
                    info = []
                    data = {}  # to hold axis#: list of numbers k: v pairs
                elif line.startswith('atom'):
                    # Warning: except for the 'name' statement, which must be
                    # at the beginning, these lines are the only ones where the
                    # order is important.  The first 'atom' line has the
                    # selection of the first atom in the first angle, and so on.
                    selstrings.append('('+line.split(None, 1)[1]+')') 
                elif line.startswith('info'):
                    info.append(line.split(None, 1)[1])
                elif line.startswith(('axis', 'energy')):
                    # Allows for the splitting of a number line into several 
                    # with the same statement (e.g., 'axis1' or 'energy').
                    words = line.split()
                    if 'axis' in line: axis = int(words[0][4:])
                    if 'energy' in line: axis = 999  # the largest value
                    data.setdefault(axis, []).extend([float(s) for s in
                                                      words[1:]])
                else:
                    msg = 'Unrecognized statement in database file %s'% filename
                    raise Exception(msg)
        # Set up the last term in filename.
        pot = setup_pot(name, '\n'.join(info), data, selstrings, verbose)
        pots.append(pot)
        potname = pot.potName()
        potnames[potname] = potnames.setdefault(potname, 0) + 1
        result.append(pots) # append terms to returned potlist

    # Write counts.
    for (potname, count) in sorted(potnames.items()):
        print('Number of %s terms: %i' % (potname, count))
    
    return result    
