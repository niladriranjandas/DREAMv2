"""Module to help set up statistical torsion angle, interpolated potentials.

Although the provided tools are general and can be used for setting up any such
potential (e.g. user-provided), this module is specially suited for setting up
the so-called torsionDB potential for protein, RNA, and DNA.  A particularly
useful function is create_TorsionDBPot, which generates the potential.  

When using this module, please cite:

For protein applications:

Bermejo, G.A., Clore, G.M., and Schwieters, C.D. (2012).  Smooth statistical
torsion angle potential derived from a large conformational database via
adaptive kernel density estimation improves the quality of NMR protein
structures.  Protein Sci. 21, 1824-1836.

For RNA applications:

Bermejo, G.A., Clore, G.M., and Schwieters, C.D. (2016).  Improving NMR
structures of RNA.  Structure 24, 806-815.

"""

import os
import math  # used by function find_minima
import glob
import sys

import torsionInterpolPot
import potList
import cdsVector
import simulation
import atom 
import atomSel
import atomSelLang # used by create_TorsionDBPot
import simulationTools
import xplorSimulation

import protocol           # used by function find_minima 
import psfGen             # used by function find_minima 
import ivm                # used by function find_minima
import torsionDBTools     # used by function find_minima
import densityEstimation  # used by function find_minima
import torsionTools       # used by function find_minima
import dihedral           # used by function find_minima
import selectTools        # used by function find_minima and setup_pot
import torsionInterpolPotTools  # used by function setup_pot
import repelPotTools  # for create_RepelPot14
import xplor




# Location of database directories within Xplor-NIH directory.
# (Added TORSIONS environment variable in ../bin/xplor.in)
databases = {'protein': os.environ['TORSIONS']+'/protein',
             'rna': os.environ['TORSIONS']+'/rna',
             'dna': os.environ['TORSIONS']+'/dna'}


# List to be filled with torsions in the system under consideration.
# A torsion is here denoted by a tuple (idi, idj), where idi and idj are the
# indices of two atoms bonded to each other via the bond around which the
# torsion is defined.  Torsions are added by get_torsions() and removed by
# setup_pot().
torsions = []


# Dict to be filled with idi:[idj, ..., idn] key:value pairs, where idx is the
# index of atom x, and the listed indices represent atoms covalently bound to
# idi (i.e., idi's "neighbors").  All atoms in the system under consideration
# appear as keys.  The dict is populated by get_torsions() and remains unchanged
# after that.
neighbors = {}


def getsystems():
    return list(databases.keys())

def getdblocation(system):
    return databases.get(system)



def getdbfilepath(system, database='default', ext='.dat'):
    """Return a list with the path(s) to database file(s).

    database is a string with the name of either a database file or a directory
    containing one or more database files.  Alternatively, a sequence with
    several such names can be provided as the database argument.  The system is
    that to which the database applies to (e.g., 'protein'), and is specified by
    the system argument; its possible values can be found as follows:

        import torsionDBPotTools  
        print( torsionDBPotTools.getsystems() )  # prints valid systems

    If database is not provided, a list with the path(s) of default database
    file(s) for the system is returned.

    The path of a name (file or directory name) in database is established by
    first searching for the name in the working directory, and then in a system-
    specific location within the Xplor-NIH directory; such location can be found
    as follows:

        print( torsionDBPotTools.getdblocation('protein') ) # prints location
                                                      # of protein databases in
                                                      # Xplor-NIH directory

    If the searched name is that of a file, its path is appended to the returned
    list.  If the searched name is that of a directory, the paths of all files
    within such directory are appended to the returned list, as long as their
    file extensions are that of the ext argument.

    """
    if not getdblocation(system):
        raise Exception("system '%s' is unsupported" % system)

    if type(database) is str: database = [database]
    result = []
    for name in database:
        temp = []
        # First, search in current directory.
        if os.path.isfile(name): temp.append(name)
        if os.path.isdir(name): temp.extend(glob.glob(name+'/*'+ext))
        # If nothing found, then search in Xplor-NIH directory.
        if not temp:
            longname = getdblocation(system) + '/' + name
            if os.path.isfile(longname): temp.append(longname)
            if os.path.isdir(longname): temp.extend(glob.glob(longname+
                                                              '/*'+ext))
        # If nothing found raise error, otherwise add findings to result.
        if not temp:
            msg = "No files associated to name '%s' in current dir or in %s" \
            % (name, getdblocation(system))
            raise Exception(msg)
        else:
            result.extend(temp)

    return result
    

def del_resid_from_metasel(metaselection, word=''):
    """Return a proper XPLOR selection string.

    Remove the 'resid word' specification from the input metaselection so that
    the remaining, returned string is a proper XPLOR selection string.

    """
    if not word: raise Exception("Specify word that follows 'resid'")
        
    # Convert any internal multiple-character spacing to single spacing
    # (needed for successful string search below).
    metaselection = ' '.join(metaselection.split())
    
    # We assume the following 3 ways in which 'resid word' may be embedded
    # in the metaselection string:
    #  ... and resid word           (at the end)
    #  ... and resid word and ...   (in the middle)
    #  resid word and ...           (at the beginning)

    # This handles the end and middle cases.
    result = metaselection.replace('and resid %s' % word, '')

    # In the case of beginning format, the above replace did nothing.
    # This replace handles the beginning case.
    return result.replace('resid %s and' % word, '')    


def metasels_central_residue_selstring(metaselections):
    """
    Find first metaselection in metaselections that is of the kind
    'sel and resid _RESID', and return sel.
    
    """
    if type(metaselections) is str: metaselections = [metaselections]

    for metasel in metaselections:
        # Avoid metaselections that don't select atoms in the residue of
        # interest (resid _RESID; the "central" residue).
        if not ('_RESID-' in metasel or '_RESID+' in metasel):
            return del_resid_from_metasel(metasel, '_RESID')


def metasels_atoms_offsets(metaselections, sim):
    """
    For each metaselection in metaselections, generically written
    'sel and resid _RESID+/-X', return (<m atomSel>.AtomSel(sel, sim), +/-X).

    """
    if type(metaselections) is str: metaselections = [metaselections]

    result = []
    for metasel in metaselections:

        # Find offset.
        words = metasel.split('_RESID')
        
        if len(words) != 2:
            msg = 'There must be exactly one _RESID word in metaselection: %s'
            raise Exception(msg % metasel)

        try:
            offset = int(words[1].split()[0])
        except ValueError:                     
            offset = 0  # _RESID not followed by +/-integer
        
        # Create atomSel.AtomSel instances.
        if offset > 0:
            selstring = del_resid_from_metasel(metasel, '_RESID+%i' % offset)
        elif offset < 0:
            selstring = del_resid_from_metasel(metasel, '_RESID%i' % offset)
        else:
            selstring = del_resid_from_metasel(metasel, '_RESID')

        # Prepare list to return.
        result.append((atomSel.AtomSel(selstring, sim), offset))

    return result

    
#@profile
def setup_pot(name, data, selection, verbose):
    """Return <m potList>.PotList with instances of  N-dimensional
    <m torsionInterpolPot>.TorsionInterpolPotND potential terms.

    The name argument (a string) is the name given to the returned
    <m potList>.PotList.  data is a sequence of PotData instances with the
    information required to create and set up the terms.  The selection argument
    is an <m atomSel>.AtomSel instance specifying the atoms subjected to the
    terms.  If verbose is True, some information is printed about the setup
    process.

    As a side effect, each torsion angle affected by the terms is removed from
    the global torsions list.

    """
    # Simulation.
    sim = selection.simulation()

    # A metaselection associated with a potential term has the following general
    # format:
    #
    # 'sel and resid _RESID+/-X'
    #
    # where sel is a standard XPLOR selection string that excludes any segid or
    # resid specifications, and involves one or more atom names.  X is an
    # integer that indicates a sequence offset relative to the "central" residue
    # (for which X=0).

    # Generate list with (segid, [(resid, resname), ...]) for residues that
    # match the input selection and the central residue type associated with the
    # terms being considered.

    # Here, we assume that all terms share the same central residue type, so
    # the metaselections associated with any term will do.
    selstring = metasels_central_residue_selstring(data[0].metaselections)

    # Avoid generating time-consuming atomSel.AtomSel instances.
    if selection.string() != 'all':  # different from default
        atomsel = atomSel.AtomSel(selstring, sim)
        intersection = atomSel.intersection(atomsel, selection)
        segments = list(selectTools.getSegsResidues(intersection, sim).items())
    else:
        segments = list(selectTools.getSegsResidues(selstring, sim).items())
        pass

    # Generate terms (applied to nothing yet, i.e., no restraints).

    if type(data) not in (list, tuple): data = [data]

    pots = TDBPotl(name)  # to hold terms 
    pots.resetPotName('TDBPotl')
    
    atoms_offsets = []
    
    for datum in data:
        
        # Potential term (applied to nothing yet, i.e., no restraints).
        pot = torsionInterpolPotTools.setup_pot(datum.name, datum.info,
                                                datum.values, sim=sim)
        # Set the term's threshold.
        if datum.elevel is not None: pot.setThreshold(datum.elevel)

        # Append to pots.
        pots.append(pot)

        # Given [metaselection1, ...] <-> ['sel1 and resid _RESID+/-X1', ...]
        # for current term, create [(atomSel.AtomSel(sel1, sim), +/-X1)), ...]
        # and append to atoms_offsets list.
        atoms_offsets.append(metasels_atoms_offsets(datum.metaselections, sim))
        
    # Generate restraint selection strings and add restraints to the terms;
    # remove associated torsions from global torsions list.

    # Use of simple, single-atom selections (e.g., 'ATOM "   A" 3 C'), generated
    # by the intersection strategy in the inner loop below, is significantly
    # faster than using "full" selections (e.g., 'segid "   A" and resid 3 and
    # name C').

    # Loop thru central residues.
    for segment in segments:
        
        segid = segment[0]
        segidSel = atomSel.AtomSel('segid "%s"' % segid, sim)
            
        for (resid, resname) in segment[1]:

            residSel = atomSel.AtomSel('resid %i' % resid, sim)
            
            # Generate restraint name.
            if segid.strip(): # assign to variable in outer loop to save time?
                restname = '%s:%s%i' % (segid.strip(), resname, resid) 
            else:
                restname = '%s%i' % (resname, resid)
            
            success = True  # generated valid selection string?

            all_selstrings = [] # for selection strings from all terms

            # Loop thru terms.
            for idx in range(len(pots)):

                selstrings = [] # for single-atom selection strings defining 
                                # the torsion(s) in current term
            
                # (Effectively) Loop thru current term's metaselections.
                # Try generating valid atom selection strings.
                for (atoms, offset) in atoms_offsets[idx]:
                   
                    atom = atomSel.intersection(atoms, segidSel,
                                                createString=False)
                    if offset:
                        atom = atomSel.intersection(atom, 'resid %i' %
                                                    (resid+offset),
                                                    createString=False)
                    else:
                        atom = atomSel.intersection(atom, residSel,
                                                    createString=False)                        

                    # atom above must represent a single atom.

                    if len(atom) > 1:  # for some reason (?)
                        msg = 'Selection of one atom required for ' + \
                              'potential:restraint %s:%s - %i attempted'
                        raise Exception( msg % (pots[idx].instanceName(),
                                                restname, len(atom)) )
                    
                    # Atom may not exist (e.g., when trying to define phi at
                    # the N-terminal residue of a protein).
                    try:
                        selstrings.append('(ATOM "%s" %i %s)' % (segid,
                                                                 resid+offset,
                                                            atom[0].atomName()))
                    except IndexError:  # when len(atom) = 0
                        if verbose:
                            # If invalid selection string for one term, no term
                            # will affect current central residue (see below).
                            msg = '%s NOT applied to %s'
                            for x in range(len(pots)):
                                print(msg % (pots[x].instanceName(), restname))
                        success = False
                        break # from metaselection loop

                # If invalid selection string for current term, don't bother
                # with other terms: Consider next central residue.
                if not success: break  # from term loop
                    
                all_selstrings.append(selstrings)  # gather selection strings

            # If all selection strings for all terms are valid:
            if success:

                for idx in range(len(pots)):

                    selstrings = all_selstrings[idx]
                    
                    # Add restraint.
                    restraint = 'name %s\nassign %s' % (restname,
                                                        ' '.join(selstrings))
                    pots[idx].addRestraints(restraint)                                                 
                    if verbose:
                        msg = '%s applied to %s'
                        print(msg % (pots[idx].instanceName(), restname))

                    # Remove associated torsions from global torsions list.
                    for i in range(len(selstrings)//4):
                        four_selstrings = selstrings[i*4: i*4+4]
                        a = atomSel.AtomSel(four_selstrings[1])[0].index()
                        b = atomSel.AtomSel(four_selstrings[2])[0].index()
                        # A torsion may be involved in more than one term;
                        # trying to remove it from torsions after it's gone 
                        # will yield ValueError.
                        try:
                            if a < b:
                                torsions.remove((a, b))
                            else:
                                torsions.remove((b, a))
                        except ValueError:
                            pass
    return pots


def get_impropers(psfstring):
    """Return a list with lists of impropers from the input PSF file.

    The PSF file is input as the string psfstring.  Each improper in the
    returned list is a list with the four atom IDs which define the improper.
    Each ID is a string with the ID number given in the PSF file.

    """
    (head, sep, tail) = psfstring.partition('!NIMPHI: impropers')

    result = []

    tail = tail.lstrip()  # remove leading spaces
    lines = tail.splitlines(True)
    for line in lines:
        if not line.isspace():
            words = line.split()
            (first, second) = (words[:4], words[4:])
            result.append(first)
            result.append(second)
        else:
            break  # reached end of impropers section
    return result


def define_cispro_impropers(psfstring, prePRO_residues):
    """Return a list with improper definitions for each PRO assuming
    cis peptide. 

    The specific molecular system is specified via a string with the
    contents of a PSF file (psfstring).  prePRO_residues is a list of
    (segment id, residue number) tuples (string, int types) indicating
    the preproline residues in the system.

    For each proline in the input system a list is made of the form (all such
    lists are packed into a list and returned):

    [(segid, resid), improper1, improper2, improper3]

    where segid (string) and resid (int) specify the segment ID and the residue
    number, respectively, and improperk (k = 1, 2, 3) defines a unique improper
    that the residue would have if it were in the cis conformation.  The
    definition of impropers is done in terms of atom IDs (see docstring of
    get_impropers function for more details).  For example,

    improper1 = [-C, -CA, +N, +CA]

    where "+" indicates the proline and "-" the residue preceeding it (e.g.,
    "+CA" is the ID of the CA of the proline.  The impropers are those "added"
    by the CIPP patch in file protein.top.
    
    """
    minus = {}
    plus = {}
    for (segid, resid) in prePRO_residues:
        minus[(segid, resid)] = {}
        plus[(segid, resid+1)] = {}

    # Read 'NATOM' PSF section.
    (head, sep, tail) = psfstring.partition('!NATOM')

    tail = tail.lstrip()  # remove leading spaces
    lines = tail.splitlines(True)
    for line in lines:
        if not line.isspace():
            words = line.split()
            index = words[0]
            try:  # empty segid field 
                resid = int(words[1])  # will fail if segid is words[1]
                segid = ''
                name = words[3]
            except ValueError:  # segid provided
                segid = words[1]
                resid = int(words[2])
                name = words[4]

            if (segid, resid) in list(minus.keys()):
                if name in ('CA', 'C', 'O'):
                    minus[(segid, resid)][name] = index

            if (segid, resid) in list(plus.keys()):
                 if name in ('N', 'CA', 'CD'):
                    plus[(segid, resid)][name] = index
        else:
            break  # reached end of NATOM section

    # Unique impropers in cis PRO (snippets from CIPP patch):
    # add improper -C -CA +N +CA
    # add improper -C -O  +N  +CA
    # add improper +CA +N  +CD -C

    result = []
    for (segid, resid) in list(minus.keys()):

        try:
            improper1 = [minus[(segid, resid)]['C'],
                         minus[(segid, resid)]['CA'],
                         plus[(segid, resid+1)]['N'],
                         plus[(segid, resid+1)]['CA']]

            improper2 = [minus[(segid, resid)]['C'],
                         minus[(segid, resid)]['O'],
                         plus[(segid, resid+1)]['N'],
                         plus[(segid, resid+1)]['CA']]
            
            improper3 = [plus[(segid, resid+1)]['CA'],
                         plus[(segid, resid+1)]['N'],
                         plus[(segid, resid+1)]['CD'],
                         minus[(segid, resid)]['C']]

            result.append([(segid, resid+1), improper1, improper2, improper3])
        except KeyError as msg:
            print('define_cispro_impropers: Warning: missing atom in', end=' ')
            print('segid "%s", resid %d:' % (segid,resid+1), end=' ')
            print(msg)
            pass

    return result



def get_cispro_residues(prePRO_residues, simulation):
    """Return a list with (segid, resid) for cis prolines.

    segid and resid are the segment id (string) and the residue number (int),
    respectively, of a cis proline in the currently loaded coordinates.
    Input prePRO_residues is a list similar to the returned one, except that it
    corresponds to the preproline residues.  If no cis prolines are found an
    empty list is returned.

    """
    # If no preprolines (i.e., no cis prolines) return an empty list.
    # This is the desired behavior.  It also has a practical reason: sometimes
    # we don't have a psf file, for example, when generating an energy surface.
    # Without this step, this function will write the psf anyway (an empty one)
    # and define_cispro_impropers() will try to read it.  
    if not prePRO_residues: return []
    
    # We want xSim in order to use its fastCommand method below.
    # fastCommand is fast because it doesn't link the topology info in the
    # old XPLOR interface and the C++ side.  For that reason fastCommand can be
    # dangerous but since we are not changing anything here, it's safe to use.
    xSim = xplorSimulation.getXplorSimulation(simulation)

    outputState = xSim.disableOutput()  # quiet XPLOR printout (old interface)
    
    # Write temporary PSF file.
    psffilename = simulationTools.mktemp('psf')
    xSim.fastCommand("write psf output = %s end" % psffilename)

    xSim.enableOutput(outputState)  # restores verbosity of old interface

    # Read PSF file.
    infile = open(psffilename, 'r')
    psfstring = infile.read()
    infile.close()

    # Delete temporary PSF file.
    os.remove(psffilename)

    # Get definitions of all impropers in terms of atom IDs.
    impropers = get_impropers(psfstring)

    # Define impropers that would be associated with cis PRO.
    cis_impropers = define_cispro_impropers(psfstring, prePRO_residues)

    result = []
    for item in cis_impropers:
        putatives = item[1:]  # the first offset contains (segid, resid)
        count = 0
        for putative in putatives:
            reverse = putative[:]
            reverse.reverse()  # either ordering in the PSF is equivalent
            if putative in impropers or reverse in impropers: count += 1
        if count == 3: result.append(item[0]) # all 3 impropers found =>
                                              # current PRO is cis
    return result
                


def get_torsions(sim):
    """Populate the global torsions list with all torsion angles in the PSF.

    A torsion angle is denoted by a tuple (idy, idz), where idy and idz are the
    indices of two atoms bonded to each other via the bond around which the
    torsion is defined.

    As a side effect, the global neighbors dictionary is populated, which
    contains idi:[idj, ..., idn] key:value pairs, where idx is the index of atom
    x, and the listed indices represent atoms covalently bound to idi (i.e.,
    idi's "neighbors").  All atoms in the PSF appear as keys.

    """
    for idx in range(sim.numBonds()):
        (a, b) = sim.bondPairByID(idx)
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)

    for idx in range(sim.numBonds()):
        (a, b) = sim.bondPairByID(idx)
        if len(neighbors[a]) > 1 and len(neighbors[b]) > 1:
            if a < b:   # to facilitate subsequent search
                torsions.append((a, b))
            else:
                torsions.append((b, a))


class TDBPotl(potList.realPotList):  # potList.PotList yields error message
    """A special potList.Potlist to hold torsion angle database potential terms.

    Methods defined here are intended for violation analysis.
    
    """
    def __init__(self, name, noPotSummary=False): 
        potList.realPotList.__init__(self, name)        
        self.noPotSummary = noPotSummary    
        
    def rms(self):
        """Not implemented."""
        return -1

    def numRestraints(self):
        """Return the number of restraints of one term in the sequence."""
        return self[0].numRestraints()

    def violations(self):
        """Return number of residues outside favored torsion angle regions."""
        # Restraint names set by setup_pot depend only on the residue, not on
        # the specific term; we take advantage of this here.
        names = []
        for term in self:
            for restraint in term.restraints():
                if restraint.violated(): names.append(restraint.name())
        return len(set(names))



class PotData():
    """Simple data structure to collect data that will be associated with an
    N-dimensional <m torsionInterpolPot>.TorsionInterpolPotND potential term.

    Attribute description:

    name (a string) is a name that may be given to the term.
    
    metaselections is a list of strings, each a metaselection for an atom that,
    along with the rest, define the torsion angle(s) applicable to the term.

    info is a comment string that can be retrieved from the term via the
    comment method.

    values is a dictionary that contains all the numerical data of the term.
    The keys are integers: for example, 1 labels the axis that corresponds to
    the first angle defined in the metaselections list.  The largest key is used
    to identify the energy values.  The numbers associated with the axes and
    energy are provided by the dict values in the form of a list.

    elevel specifies an energy level above which the corresponding torsion
    angle(s) in the structure under study will be flagged as "violated" (note
    that any scaling on the potential will not affect this functionality;
    unscaled energies are considered).     

    """
    def __init__(self, name, metaselections, info, elevel, values):
        self.name = name
        self.metaselections = metaselections
        self.info = info
        self.elevel = elevel
        self.values = values
        


#@profile
from trace import notrace_decorate
@notrace_decorate
def create_TorsionDBPot(name='', selection='all', system='protein',
                        database='default', threshold='99.95', ext='.dat', 
                        verbose=False):
    """Return a <m potList>.PotList of statistical torsion angle potential terms.

    Sets up torsion angle potential terms of the type defined in the
    <m torsionInterpolPot> module, based on statistics collected from a
    database.  The name argument gives the name to the returned
    <m potList>.PotList.  The type of system under consideration (e.g.,
    'protein', 'rna', 'dna') is specified by the system argument (a string).
    (A special value for system, 'repel14', not included above is discussed at
    the end.)
    
    The database argument is a string with the name of either a database file or
    a directory containing one or more database files (in the latter case, only
    files with the extension given by the ext argument are considered).
    Alternatively, a sequence with several such filenames can be provided as the
    database argument.  If no database files are specified, the default files
    for the corresponding system are used.  The path of a database file or
    directory is established (with getdbfilepath function) by first searching
    in the working directory, and then in a system-specific location within the
    Xplor-NIH directory; such location can be found as follows:

        print( torsionDBPotTools.getdblocation('protein') ) # prints location
                                                      # of protein databases in
                                                      # Xplor-NIH directory

    For documentation on the format of database files see the README.format
    file in the databases/torsions directory within the Xplor-NIH directory.
    The selection argument specifies the atoms to subject to the database
    potential (all atoms by default); it can be an XPLOR-style selection string
    or an <m atomSel>.AtomSel instance.  Note that certain torsion angles, such
    as protein phi (and psi), involve atoms from more than one residue.  Thus,
    for example, selection='resname PRO' would not enforce a hypothetical 1D PRO
    phi potential on prolines not immediately preceeded by a proline because the
    C atom from the pre-proline residue, needed to define phi, is not selected.
    
    The threshold argument is a string (possible values: '98' or '99.95') that
    specifies a percentage used to flag torsion angle combinations in the
    structure under study whose associated energy is above that of the given
    percentage of residues in the database.

    If verbose is True, some information is printed about the setup process.

    Each database file (explicitly or implicitly specified by the database
    argument) may specify one or more potential terms, which are packed into a
    <m potList>.PotList.  The latter is appended to the returned potList, which
    thus consists of a potList of potLists.

    ----------------------------------------------------------------------------

    In addition to specifying the system under study (e.g., 'protein', 'rna',
    'dna'), which configures a statistical torsion angle potential of the
    appropriate kind, the system argument may be set to 'repel14'.  In this
    special case, the function returns a <m repelPot>.RepelPot instance (i.e.,
    a repulsive-only potential) involving atoms of the type of A and D:

        A---B---C---D    
    
    where the torsion angle defined around the B---C bond is not already
    affected by a statistical torsion angle potential returned by this function
    (i.e., by using 'protein', 'rna', or 'dna' system values).

    """    
    if system != 'repel14' and system not in getsystems():
        msg = "Wrong system name provided ({})\n".format(system)
        msg += "\tchoose from one of these: {}".format(', '.join(getsystems())+
                                                       ",'repel14'")
        raise Exception(msg)
    
    # Convert input selection to an atomSel.AtomSel.
    selection = selectTools.convertToAtomSel(selection)

    # Simulation.
    sim = selection.simulation()
 
    # If the (global) torsions list is empty...
    if not torsions:
        get_torsions(sim) # ... populate it with all torsion angles in the PSF

    if system == 'repel14':
        return create_RepelPot14(sim, verbose=verbose)
        
    # Protein-specific definitions.
    if system == 'protein':
        # Define atom selection name "prePRO".  
        prePRO = atomSel.AtomSel('name C and bondedto (name N and resname PRO)',
                                 sim)
        prePRO_residues = [(x.segmentName(), x.residueNum()) for x in prePRO]
        selstring = ['(segid "%s" and resid %i)' % (x, y) for (x, y) in
                     prePRO_residues]
        prePRO = atomSel.AtomSel(' or '.join(selstring), sim)
        atomSelLang.setRecallSelection('prePRO', prePRO)

        # Define atom selection name "cisPRO".  
        cisPRO_residues = get_cispro_residues(prePRO_residues, sim)
        selstring = ['(segid "%s" and resid %i)' % (x, y) for (x, y) in
                     cisPRO_residues]
        cisPRO = atomSel.AtomSel(' or '.join(selstring), sim)
        atomSelLang.setRecallSelection('cisPRO', cisPRO)

    # Set up potential.
    result = potList.PotList(name)
    filenames = getdbfilepath(system, database, ext)
    for filename in filenames:

        if verbose: print('Processing', filename)

        # Create a potlist for all terms within the current file.
        potlistname = filename.split('/')[-1] # get filename without path
        potlistname = potlistname.split('.')[0]  # discard filename extension

        # To pack all data associated with the different terms in current file.
        data = []
        
        infile = open(filename, 'r')
        contents = infile.read()
        infile.close()

        name = ''  # now the name of the potential term (not the result potlist)
        lines = contents.splitlines()
        foundThreshold = False
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#'):
                if line.startswith('name'):
                    if name: # the start of a new term; set up the previous one
                        data.append(PotData(name, metaselections,
                                            '\n'.join(info), elevel, values))
                    name = line.split()[1]
                    metaselections = []
                    info = []
                    elevel = None 
                    values = {}  # to hold axis no.:list of numbers k:v pairs
                elif line.startswith('atom'):
                    # Warning: except for the 'name' keyword, which must be
                    # at the beginning, these lines are the only ones where the
                    # order is important.  The first 'atom' line has the meta-
                    # selection of the first atom in the first angle, and so on.
                    # What follows the 'atom' keyword is called a metaselection
                    # because it's not necessarily a proper XPLOR selection due
                    # to, e.g., the '_RESID' word.
                    metaselections.append(line.split(None, 1)[1])
                elif line.startswith('info'):
                    info.append(line.split(None, 1)[1])
                elif line.startswith('elevel'): 
                    words = line.split()                
                    if words[1] == threshold:
                        elevel = float(words[2])
                        foundThreshold=True
                        pass
                    pass
                elif line.startswith('axis'):
                    # Allows for the splitting of the line into several starting
                    # with the same keyword.
                    words = line.split()
                    axis = int(words[0][4:])
                    values.setdefault(axis, []).extend(words[1:])
                elif line.startswith('energy'):
                    # Allows for the splitting of the line into several starting
                    # with the same keyword.
                    words = line.split()
                    values.setdefault(999, []).extend(words[1:])
                else:
                    msg = 'Unrecognized keyword in database file %s' % filename
                    raise Exception(msg)
                pass
            pass
        
        # At this point we've gone thru all terms in the current file.

        # At least one term in the current file must have the specified
        # threshold (and associated elevel) defined.  (The threshold isn't
        # necessarily applicable to all terms in the file.)
        if threshold != None and not foundThreshold:
            raise Exception('threshold value %s not found in %s' % (threshold,
                                                                    filename))
        # Pack data associated with the last term in current file.
        data.append(PotData(name, metaselections, '\n'.join(info), elevel,
                            values))

        # Set up terms and append to returned potlist.
        result.append(setup_pot(potlistname, data, selection, verbose))
        pass

    # Write counts.
    from simulationWorld import SimulationWorld_world
    if SimulationWorld_world().logLevel() != 'none' or verbose:
        names = [x.instanceName() for x in result]
        names.sort()
        print('Number of residues affected by each database file:')
        for name in names:
            print('%s: %i' % (name, result[name][0].numRestraints()))
        
    return result


##def create_RepelPot14(sim, verbose=False):
##    """Return a <m repelPot>.RepelPot instance.
##
##    The returned potential contains 1-4 interactions between atoms bound to
##    idi and idj for each (idi, idj) tuple in global torsions list.
##
##    """
##    # May want to add selection to exclude from consideration (e.g., PHE rings)
##    # to speed up things.
##    
##    if verbose:
##        print 'torsionDBPotTools adding 1-4 repulsive interactions between:'
##    
##    pairs = []
##    
##    for (a, b) in torsions:
##
##        # Atoms bound to a, excluding b.
##        a_side = neighbors[a][:]  # copy (don't modify global neighbors dict) 
##        a_side.remove(b)
##
##        # Atoms bound to b, excluding a.
##        b_side = neighbors[b][:]  # copy (don't modify global neighbors dict)
##        b_side.remove(a)
##        
##        if a_side and b_side:
##            for x in a_side:
##                pair = (atomSel.AtomSel(x, sim), atomSel.AtomSel(b_side, sim))
##                pairs.append(pair)
##                if verbose:
##                    print '(%s) and: %s' % (pair[0][0].string(),
##                            ', '.join(['('+x.string()+')' for x in pair[1]]))
##
##    # Return repelPot.RepelPot instance.
##    # repel argument not specified: it's automatically set according to the
##    # parameters used.
##    return repelPotTools.create_RepelPot('repel14', use14=True, selPairs=pairs)


def create_RepelPot14(sim, verbose=False):
    """Return a <m repelPot>.RepelPot instance.

    The returned potential contains 1-4 interactions between atoms bound to
    idi and idj for each (idi, idj) tuple in global torsions list.

    """
    # May want to add selection to exclude from consideration (e.g., PHE rings)
    # to speed up things.
    
    if verbose:
        print('torsionDBPotTools adding 1-4 repulsive interactions between:')
    
    pairs = []
    
    for (a, b) in torsions:

        # Atoms bound to a, excluding b.
        a_side = neighbors[a][:]  # copy (don't modify global neighbors dict) 
        a_side.remove(b)

        # Atoms bound to b, excluding a.
        b_side = neighbors[b][:]  # copy (don't modify global neighbors dict)
        b_side.remove(a)
        
        for x in a_side:
            for y in b_side:
                # To avoid redundant interactions in rings.
                if (x, y) not in pairs and (y, x) not in pairs:
                    pairs.append((x, y))
                    
    pairs = [(atomSel.AtomSel(x, sim), atomSel.AtomSel(y, sim)) for (x, y) in
             pairs]

    if verbose:
        for pair in pairs:
            print('(%s) and (%s)' % (pair[0][0].string(), pair[1][0].string()))

    # Return repelPot.RepelPot instance.
    # repel argument not specified: it's automatically set according to the
    # parameters used.
    return repelPotTools.create_RepelPot('repel14', use14=True, selPairs=pairs)




def find_minima(sequence, resid=2, cispeptide=[], database=[], system='protein',
                exclude='omega', npoint=30, tolerance=1.0, flag=0, nmin=3,
                numSteps=500, dEPred=0.001, printInterval=0):
    """Return a list with (energy, minimum) tuples sorted by energy value.

    Each tuple in the returned list contains the energy value of a local minimum
    and a <m cdsVector>.CDSVector_double with the corresponding coordinates in
    torsion angle space.  The sequence argument is a string with the residue
    sequence (e.g., 'GLY ALA GLY'), and resid is a residue number in it.  Cis
    peptide bonds can be specified by via cispeptide, a sequence of one or more
    residue numbers (int), each indicating the residue N-terminal to the cis
    peptide bond; here a blank segment ID (segid) is assumed (if non-blank
    segids are required, the elements of cispeptide may consist of (resiue#,
    segid) pairs, segid a string).

    The minima are found by minimization of all torsion angles in resid,
    subjected to the torsion angle database potential specified by the database
    and system arguments (see create_TorsionDBPot function for their
    description).  exclude is either a string with a torsion name or a list of
    such names (for all applicable names see <m torsionDBTools> module); the
    associated torsion angles are assumed to have no effect on the overall
    energy, and are thus disregarded.  Several minimizations are performed,
    starting from torsion coordinates in a grid; the number of points in each
    dimension are specified by the npoint argument.  Two multidimensional
    torsion angle points, minimum1 and minimum2, are considered to represent
    different minima if the short angle between them is greater than tolerance;
    if such angle is between tolerance and the flag argument, one of the torsion
    angle points is flagged in the returned list as follows:

    (energy1, minimum1, [minimum2, small-angle]).

    The above is only done if flag!=0.  nmin is the number of consecutive
    minimizations to perform from a given starting point; the rest of the
    arguments are those described in function <m protocol>.initMinimize.

    """
    # Generate molecule.
    
    protocol.initStruct()  # clear previous structure info
    psfGen.seqToPSF(sequence)  # parameters loaded as side effect

    for item in cispeptide:  # handle cis peptide bonds if any
        if type(item) == int:  # only a resid is provided
            psfGen.cisPeptide(item)
        else:  # assume item is a sequence
            psfGen.cisPeptide(item[0], item[1])  # resid and segid provided

    protocol.genExtendedStructure()

    # Get atom selections associated to all torsions in specified residue.
    resname = sequence.split()[resid-1]
    torsions = torsionDBTools.get_torsion_names(resname) # list of torsion names

    # Exclude selections of torsions indicated in the exclude argument.
    if type(exclude) == str: exclude = [exclude]
    for torsion in exclude:
        torsions.remove(torsion.strip().lower())
    
    # Tuple of tuples, each with 4 atom selections defining torsions.
    selections = torsionDBTools.get_torsion_selections(resname, torsions, resid)

    # Generate a selection by merging all in selections; it will define atoms
    # to consider for the potential setup.
    selection = [' or '.join(x) for x in selections]
    selection = ' or '.join(selection)

    # Set up potential.
    potlist = create_TorsionDBPot(database=database, system=system,
                                  selection=selection)
    # IVM object for minimization.
    mint = ivm.IVM()
    protocol.torsionTopology(mint)
    protocol.initMinimize(ivm=mint, potList=potlist, numSteps=numSteps,
                          dEPred=dEPred, printInterval=printInterval)

    # Generate grid of torsion values to use as starting points in minimization.
    ndim = len(torsions) # dimensionality of torsion-angle space to be probed
    ranges = torsionDBTools.get_torsion_range(resname, torsions)
    domains = [max(x)-min(x) for x in ranges]
    points = ndim * [npoint]
    grid = densityEstimation.grid(points, ranges)

    # To get the minimized torsion values.
    dihedrals = [dihedral.Dihedral(x[0], x[1], x[2], x[3]) for x in selections]

    # Loop through each point in grid and perform minimization.
    result = []
    for point in range(grid.cols()):
        init_torsions = []  # to hold ((sel1,...,sel4), torsion_value) tuples
        for dim in range(ndim):
            init_torsions.append((selections[dim], grid[dim, point]))

        # Set torsions to initial values and minimize.
        torsionTools.setTorsions(angles=init_torsions, ivm=mint)
        
        # More than one minimization may be required (sometimes it gets stuck).
        for i in range(nmin):
            mint.run()

            # For testing.
##            import simulationTools
##            simulationTools.testGradient(potlist)

        selectTools.correctSymmetricSidechains(sel=selection)

        # Get minimized torsion values.
        newmin = [math.degrees(x.value()) for x in dihedrals]

        # Include them in the result only if not already there.
        flagged = []
        for item in result:
            oldmin = item[1]
            distance = torsionDBTools.get_smallangle(oldmin, newmin, domains)
            if distance < tolerance: break
            if distance < flag:
                flagged.append((oldmin, distance))
        else:
            if flag:
                result.append((potlist.calcEnergy(), newmin, flagged))
            else:
                result.append((potlist.calcEnergy(), newmin))

    result.sort()  # sorted by energy
    return result
    


def find_minima_parallel(sequence, resid=2, cispeptide=[], database=[],
                         system='protein', exclude='omega', npoint=30,
                         tolerance=1.0, flag=0, nmin=3, numSteps=500,
                         dEPred=0.001, printInterval=0):
    """Return a list with (energy, minimum) tuples sorted by energy value.

    (Note: this function is the same as find_minima, except that it has to be
    run in parallel.)
    
    Each tuple in the returned list contains the energy value of a local minimum
    and a <m cdsVector>.CDSVector_double with the corresponding coordinates in
    torsion angle space.  The sequence argument is a string with the residue
    sequence (e.g., 'GLY ALA GLY'), and resid is a residue number in it.  Cis
    peptide bonds can be specified by via cispeptide, a sequence of one or more
    residue numbers (int), each indicating the residue N-terminal to the cis
    peptide bond; here a blank segment ID (segid) is assumed (if non-blank
    segids are required, the elements of cispeptide may consist of (resiue#,
    segid) pairs, segid a string).

    The minima are found by minimization of all torsion angles in resid,
    subjected to the torsion angle database potential specified by the database
    and system arguments (see create_TorsionDBPot function for their
    description).  exclude is either a string with a torsion name or a list of
    such names (for all applicable names see <m torsionDBTools> module); the
    associated torsion angles are assumed to have no effect on the overall
    energy, and are thus disregarded.  Several minimizations are performed,
    starting from torsion coordinates in a grid; the number of points in each
    dimension are specified by the npoint argument.  Two multidimensional
    torsion angle points, minimum1 and minimum2, are considered to represent
    different minima if the short angle between them is greater than tolerance;
    if such angle is between tolerance and the flag argument, one of the torsion
    angle points is flagged in the returned list as follows:

    (energy1, minimum1, [minimum2, small-angle]).

    The above is only done if flag!=0.  nmin is the number of consecutive
    minimizations to perform from a given starting point; the rest of the
    arguments are those described in function <m protocol>.initMinimize.

    """
    # Generate molecule.
    
    protocol.initStruct()  # clear previous structure info
    psfGen.seqToPSF(sequence)  # parameters loaded as side effect

    for item in cispeptide:  # handle cis peptide bonds if any
        if type(item) == int:  # only a resid is provided
            psfGen.cisPeptide(item)
        else:  # assume item is a sequence
            psfGen.cisPeptide(item[0], item[1])  # resid and segid provided

    protocol.genExtendedStructure()

    # Get atom selections associated to all torsions in specified residue.
    resname = sequence.split()[resid-1]
    torsions = torsionDBTools.get_torsion_names(resname) # list of torsion names

    # Exclude selections of torsions indicated in the exclude argument.
    if type(exclude) == str: exclude = [exclude]
    for torsion in exclude:
        torsions.remove(torsion.strip().lower())
    
    # Tuple of tuples, each with 4 atom selections defining torsions.
    selections = torsionDBTools.get_torsion_selections(resname, torsions, resid)

    # Generate a selection by merging all in selections; it will define atoms
    # to consider for the potential setup.
    selection = [' or '.join(x) for x in selections]
    selection = ' or '.join(selection)

    # Set up potential.
    potlist = create_TorsionDBPot(database=database, system=system,
                                  selection=selection)
    # IVM object for minimization.
    mint = ivm.IVM()
    protocol.torsionTopology(mint)
    protocol.initMinimize(ivm=mint, potList=potlist, numSteps=numSteps,
                          dEPred=dEPred, printInterval=printInterval)

    # Generate grid of torsion values to use as starting points in minimization.
    ndim = len(torsions) # dimensionality of torsion-angle space to be probed
    ranges = torsionDBTools.get_torsion_range(resname, torsions)
    domains = [max(x)-min(x) for x in ranges]
    points = ndim * [npoint]
    grid = densityEstimation.grid(points, ranges)

    # To get the minimized torsion values.
    dihedrals = [dihedral.Dihedral(x[0], x[1], x[2], x[3]) for x in selections]

    # Loop through each point in grid and perform minimization.
    # Use of the simulationTools module for structure calculation allows
    # parallelization, which might be handy in large dimensions.
    def calcOneStructure(loopInfo):

        init_torsions = []  # to hold ((sel1,...,sel4), torsion_value) tuples
        for dim in range(ndim):
            init_torsions.append((selections[dim],
                                  grid[dim, loopInfo.structureNum()]))

        # the structure number is the same as the column number in the grid. 

        # Set torsions to initial values and minimize.
        torsionTools.setTorsions(angles=init_torsions, ivm=mint)
        
        # More than one minimization may be required (sometimes it gets stuck).
        for i in range(nmin):
            mint.run()

            # For testing.
##            import simulationTools
##            simulationTools.testGradient(potlist)

        selectTools.correctSymmetricSidechains(sel=selection)
   
        # Get minimized torsion values.
        minimum = [math.degrees(x.value()) for x in dihedrals]

        # Pack it with the energy.
        global ene_minimum  # need to access ene_minimum outside function
        ene_minimum = (potlist.calcEnergy(), minimum)


    simulationTools.StructureLoop(numStructures=grid.cols(),
                                  structLoopAction=calcOneStructure).run()

    # At this point we have one ene_minimum per structure.
    ene_minima = xplor.p_comm.collect(ene_minimum) # collect into a list
    if xplor.p_processID != 0: sys.exit(0) # kill all but process 0

    # Remove repetitions and flag close minima.
    result = []
    for item1 in ene_minima:
        flagged = []
        for item2 in result:
            distance = torsionDBTools.get_smallangle(item1[1], item2[1],
                                                     domains)
            if distance < tolerance: break
            if distance < flag:
                flagged.append((item2[1], distance))
        else:
            if flag:
                result.append((item1[0], item1[1], flagged))
            else:
                result.append((item1[0], item1[1]))            
            
    result.sort()  # sorted by energy
    return result



def create_Terminal14Pot(name, repel=None, selection="not PSEUDO"):
    """Return a <m repelPot>.RepelPot instance with 1-4 interactions of terminal
    protonated groups.

    Statistical torsion angle potentials are designed to replace parametric
    terms like the <l https://nmr.cit.nih.gov/xplor-nih/xplorMan/node116.html#eq:dihedral older XPLOR dihedral potential>.
    They perform better when the van der Waals interactions between atoms
    separated by three bonds (i.e., 1-4 interactions) are turned off, as such
    interactions types are already implied in the potential.  However, by
    nature, statistical torsion angle potentials only involve torsions defined
    by heavy atoms, leaving those involving terminal protons (e.g., methyl
    groups) unrestrained.  As a result, if no 1-4 interactions are allowed for
    such groups, unrealistic eclipsed conformations may be obtained.  This
    function returns a <m repelPot>.RepelPot instance with the 1-4 nonbonded
    interactions of such terminal protonated groups.

    The argument name is the instanceName assigned to the returned
    <m repelPot>.RepelPot.  The optional repel argument (float) is the scale
    factor for the atomic radii; its optimal value depends on the version of the
    parameter file used.  The selection argument is a string or <m atomSel>.AtomSel
    object used to specify a subset of the atoms (or an alternate
    <m simulation>.Simulation.)  FIX: improve this paragraph.

    In the future this function will be deprecated in favor of a dedicated,
    parametric torsion angle potential that complements the statistical torsion
    angle term used.

    """
    selection = selectTools.convertToAtomSel(selection)
    sim = selection.simulation()

    # Selected 1-4 interactions.
    # values: [(proton, atom1, atom2, atom3), ...], where proton is, e.g., that
    # of the methyl.
    # The termini are missing!
    resTypePairs = {
        # Protein.
        'ALA' : [('HB#', 'HA', 'N', 'C')],
        'ARG' : [],
        'ASN' : [],
        'ASP' : [('HD2', 'CB', 'OD1')], # depends on protonation state
        'CYS' : [('HG', 'HB#', 'CB')],  # depends on disulfide
        'GLN' : [],
        'GLU' : [('HE2', 'CG', 'OE1')], # depends on protonation state
        'GLY' : [],
        'HIS' : [],
        'ILE' : [('HG2#', 'HB', 'CG1', 'CA'), ('HD1#', 'HG1#', 'CB')],
        'LEU' : [('HD1#', 'HG', 'CD2', 'CB'), ('HD2#', 'HG', 'CD1', 'CB')],
        'LYS' : [('HZ#', 'HE#', 'CD')],
        'MET' : [('HE#', 'CG')],
        'PHE' : [],
        'PRO' : [],
        'SER' : [('HG', 'HB#', 'CA')],
        'THR' : [('HG2#', 'HB', 'OG1', 'CA'), ('HG1', 'HB', 'CG2', 'CA')],
        'TRP' : [],
        'TYR' : [('HH', 'CE#')], # probably not important
        'VAL' : [('HG1#', 'HB', 'CG2', 'CA'), ('HG2#', 'HB', 'CG1', 'CA')],
        
        # Nucleic acids.
        # (RNA OH proton is H2', as in nucleic-3.1.)
        'GUA' : [("H2'", "H2''", "C3'", "C1'")],
        'ADE' : [("H2'", "H2''", "C3'", "C1'")],
        'PUR' : [("H2'", "H2''", "C3'", "C1'")],
        'CYT' : [("H2'", "H2''", "C3'", "C1'")],
        'THY' : [("H2'", "H2''", "C3'", "C1'"), ('H5#', 'C4', 'C6')],
        'URI' : [("H2'", "H2''", "C3'", "C1'")],
        'AP7' : [("H2'", "H2''", "C3'", "C1'")],
        'TED' : [("H2'", "H2''", "C3'", "C1'")]
                    }
    
    pairs14 = []
    import atomSel
    for atom in atomSel.AtomSel('tag',sim):
        resn = atom.residueName()
        if not resn in list(resTypePairs.keys()): continue
        resid = atom.residueNum()
        segid = atom.segmentName()
        pairs = resTypePairs[resn]
        for pair in pairs:
            a = 'segid "%s" and resid %d and (name %s)' % (segid, resid,
                                                           pair[0])
            if len(atomSel.AtomSel(a,sim)) == 0: continue
            b = 'segid "%s" and resid %d and (name %s)' % (segid, resid,
                                                    " or name ".join(pair[1:]))
            pairs14.append((a, b))

    import repelPotTools
    result = repelPotTools.create_RepelPot(name, repel=repel,
                                           use14=True, selPairs=pairs14,
                                           selection=selection)
    return result



def analyze(potlist):
    """Perform analysis of torsionDBPot terms.

    """
    potlist = simulationTools.getPotTerms(potlist, 'TDBPotl')

    if not potlist: return ''  # this goes to pdb file header

    violated = {}
    for container in potlist:  # 'container' is a potlist
        for term in container:
            for restraint in term.restraints():
                if restraint.violated():

                    # Identify the restrained residue.
                    # I could use restraint.name() to get this info, but this
                    # should be more robust.  That said, this may need to change
                    # if a new term is added that acts on more than one residue.

                    # We just look at the second atom of the first torsion
                    # affected by the term.
                    resid = restraint.angle1.atom1().residueNum()
                    resname = restraint.angle1.atom1().residueName()
                    segid = restraint.angle1.atom1().segmentName()
                    if '1D' in term.potName():
                        torsions = (restraint.angle1val,)
                    if '2D' in term.potName():
                        torsions = (restraint.angle1val, restraint.angle2val)
                    if '3D' in term.potName():
                        torsions = (restraint.angle1val, restraint.angle2val,
                                    restraint.angle3val)
                    termname = term.instanceName()
                    threshold = term.threshold()
                    energy = restraint.energy()
                    violated.setdefault(segid, []).append((resid, resname,
                                                           termname, threshold,
                                                           energy, torsions))
    segids = list(violated.keys())
    segids.sort()
    for segid in segids:
        violated[segid].sort()  # sort using residue number

    result = 'Residues that violate torsionDBPot \n\n'
    result += '%-13s %-30s %-20s %-10s %-6s \n' % ('Residue', 'Term',
                                               'Threshold energy', 'Energy',
                                               'Angles')
    for segid in segids:
        for item in violated[segid]:
            if len(item[-1]) == 1:
                angles = '(%6.3f,)' % item[-1][0]
            if len(item[-1]) == 2:
                angles = '(%6.3f, %6.3f)' % (item[-1][0], item[-1][1])
            if len(item[-1]) == 3:
                angles = '(%6.3f, %6.3f, %6.3f)' % (item[-1][0], item[-1][1],
                                                    item[-1][2])
            residue = '%i %s %s' % (item[0], item[1], segid)
            
            result += '%-13s %-30s %-20.3f %-10.3f ' % (residue, item[2],
                                                        item[3], item[4])
            result += angles + '\n'
            
    print(result)      # this goes to .viols file
    return ''  # this goes to pdb file header


simulationTools.registerTerm(analyze, 'TorsionDBPot terms', 'TDBPotl',
r"""
This term reports no additional header information.
""")
    


                    


            
            

        
    

                     
                    
