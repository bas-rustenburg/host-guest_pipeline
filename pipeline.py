# -*- coding: utf-8 -*-
"""
@author: Bas Rustenburg
"""


import sys
import argparse
import tempfile
import os
import subprocess
from glob import glob as glob
import rdkit.Chem.AllChem as Chem

# Add custom python module directory to path
sys.path.append("/home/rustenburg/Software/pythonlibs")
from mmtools.moltools import ligandtools as lt

def try_opening(parser, arg, mode="r"):
    """
    Open the file supplied to the argparser and return file object, or pass the error on to the parser.
    ARGUMENTS
    ---------
        parser (argparse.Argumentparser) - Parser that is receiving the input
        arg (string) - Filename supplied as command line input
    OPTIONAL ARGUMENTS
    ------------------
        mode (string) - Mode to be used to open the file (default: "r")

    RETURNS
    -------
        Python file object opened in the mode selected.

    REQUIREMENTS
    ------------
        argparse module (python2.7)

    TODO:
        Need a separate function for any output files?
    """
    try:
        f = open(arg, mode)
        return f
    except IOError as e:
        parser.error("Unable to open file %s: %s" % (arg, e))


def load_amber(amberhome="/home/rustenburg/Software/amber14"):
    """
    Set system variables to access Amber installation
    """

    # See if things were already set
    try:
        os.environ['AMBERHOME']
        return
    except:
        os.environ['AMBERHOME'] = amberhome

    # if AMBERHOME wasn't set, we will likely have to set PATH and
    # LD_LIBRARY_PATH as well

    # Test if PATH exists, if so append, else set it
    try:
        os.environ['PATH']  # error if non-existent
        os.environ['PATH'] = os.environ['PATH'] + ":" + amberhome + "/bin"
    except:
        os.environ['PATH'] = amberhome + "/bin"

    # Test if LD_LIBRARY_PATH exists, if so append, else set it
    try:
        os.environ['LD_LIBRARY_PATH']  # error if non-existent
        os.environ['LD_LIBRARY_PATH'] = os.environ[
            'LD_LIBRARY_PATH'] + ":" + amberhome + "/lib"
    except:
        os.environ['LD_LIBRARY_PATH'] = amberhome + "/lib"
    return


def load_oechem(oelicense="/home/rustenburg/Licenses/oe_license.txt", override=False):
    """
    Load openeye libraries and check if they're valid'
    """
    #See if import is succesful
    try:
        import openeye
        import openeye.oechem as oechem

    except ImportError,e:
        return tuple([False,e, None])

    #Make sure license can be found
    try:
        os.environ['OE_LICENSE']
        if override:
            os.environ['OE_LICENSE'] = oelicense
    except:
        os.environ['OE_LICENSE'] = oelicense

    #Check if OEChem is licensed
    if openeye.OEChemIsLicensed():
        return tuple([True, "", oechem])
    else:
        return tuple([False, "LICENSE: OEChem is unlicensed.", None])

class OpenEyePipeline(object):

    """
    Object for holding molecule data and framework for performing operations
    ARGUMENTS
    ---------
    name - str
        name of the compound
    iupac - str
        IUPAC name of the compound
    formal_charge - int
        the formal charge of the compound
    wdir - str
        working directory
    """

    def __init__(self, name="molecule", iupac="", formal_charge=0, wdir="./"):
        """
        Initialize molecular data structure

        ARGUMENTS
        ---------
            name - str
                name of the compound
            iupac - str
                IUPAC name of the compound
            formal_charge - int
                the formal charge of the compound
            wdir - str
                working directory
        """
        self.name = name
        self.iupac = iupac
        self.formal_charge = int(formal_charge)
        self.charge_assigned = False
        self.wdir = wdir
        self.mol = list()
        oechemavailable,oechemerror,self.oechem=load_oechem()

        if not oechemavailable:
            raise Exception(oechemerror)

    def create_from_iupac(self, **kwargs):
        """
        Creates a molecule from IUPAC name.
        ARGUMENTS
            **kwargs - keyword args passed directly to createMoleculeFromIUPAC.
        """
        self.mol.append(
            lt.createMoleculeFromIUPAC(
                self.iupac,
                charge=self.formal_charge,
                **kwargs))
        return self.mol

    def create_from_file(self, filename):
        """Setup a molecule from a file. """
        mol = lt.readMolecule(filename)
        self.formal_charge = lt.formalCharge(mol)
        self.iupac = lt.OECreateIUPACName(mol)
        self.mol.append(lt.expandConformations(mol, maxconfs=1))
        return self.mol

    def load_conformers(self,pattern):
        """Load molecule conformers from files"""
        self.mol = list()
        files = glob(pattern)
        for f in files:
            self.create_from_file(f)

    def gen_conformers(self,maxconfs=15):
        confs = list()
        confs.append(lt.expandConformations(self.mol[0], maxconfs=maxconfs))
        self.mol = confs
        return confs

    def write_multi(
            self,
            filename,
            outdir="./",
            ext="mol2",
            substructure_name='MOL',
            preserve_atomtypes=False):
        """Write a molecule with multiple conformers to multiple numbered files. Detects extension from automatically from ext.
        ARGUMENTS
        ---------
            filename (string) - the file to write the molecule to (type autodetected from filename)

        OPTIONAL ARGUMENTS
        ------------------
            substructure_name (String) - if a mol2 file is written, this is used for the substructure name (default: 'MOL')
            preserve_atomtypes (bool) - if True, a mol2 file will be written with atom types preserved

        RETURNS
        -------
            List of filenames

        NOTES
        -----
            Multiple conformers are written to separate files.
        """

        # Open an output stream.
        ostream = self.oechem.oemolostream()
        files = list()

        # Define internal function for writing multiple conformers to an output
        # stream.
        def _write_all_conformers(ostream, mols):
            """Write multiple conformers to single stream, but separate file"""
            for molecule in mols:
                for n, conformer in enumerate(molecule.GetConfs()):
                    filenum = "%s%s_%d.%s" % (outdir, filename, n, ext)
                    files.append(filenum)
                    ostream.open(filenum)
                    if preserve_atomtypes:
                        self.oechem.OEWriteMol2File(ostream, conformer)
                    else:
                        self.oechem.OEWriteConstMolecule(ostream, conformer)
                    ostream.close()

            return

        # Write all the conformers of a single molecule to numbered files.
        _write_all_conformers(ostream, self.mol)

        # Replace substructure name if mol2 file.
        for fn in files:
            suffix = os.path.splitext(filename)[-1]
            if (suffix == '.mol2' and substructure_name is not None):
                lt.modifySubstructureName(fn, substructure_name)

        return files



    def submit_for_pKa(self, software="MoKa", pH="7.4", **kwargs):
        """
        DUMMY FUNCTION

        PLAN
	----
            Takes an OEChem molecule and submits it to a protonation state prediction using a supplied software package
        ARGUMENTS
	---------
            molecule - an OEChem molecule object
            software - the software package desired to calculate the pKa
            pH - the desired pH for which the protonation state needs to be estimated
            kwargs - any keyworded arguments that need to be supplied to the software
        """
        pass
        return self.mol

    def charge(self, antechamber=True, **kwargs):
        """
        Assign partial charges with antechamber or openeye
        ARGUMENTS
            antechamber - Use antechamber for charges (boolean), or openeye if false.
        """
        if antechamber:
            try:
                # gives error if amberhome was not found.
                os.environ['AMBERHOME']
                for m,molecule in enumerate(self.mol):
                    chmolecule = self.assignPartialChargesWithAntechamber2(
                    molecule,
                    netcharge=self.formal_charge,
                    **kwargs)
                    self.mol[m] = chmolecule
            except KeyError:
                raise Exception(
                    "AMBERHOME was not set. Cannot detect AMBER installation: %s")

        else:
            for m,molecule in self.mol:
                self.mol[m] = lt.assignPartialCharges(molecule, **kwargs)

        self.charge_assigned = True
        return self.mol

    def assignPartialChargesWithAntechamber2(
            self,
            molecule,
            charge_model='bcc',
            judgetypes=None,
            cleanup=False,
            verbose=False,
            netcharge=None):
        """Assign partial charges to a molecule.

        ARGUMENTS
        ---------
          molecule (OEMol) - molecule for which charges are to be computed

        OPTIONAL ARGUMENTS
        ------------------
          charge_model (string) - antechamber partial charge model (default: 'bcc')
          judgetypes (integer) - if specified, this is provided as a -j argument to antechamber (default: None)
          cleanup (boolean) - clean up temporary files (default: True)
          verbose (boolean) - if True, verbose output of subprograms is displayed
          netcharge (integer) -- if given, give -nc (netcharge) option to antechamber in calculation of charges

        RETURNS
        -------
          charged_molecule (OEMol) - the charged molecule with GAFF atom types

        REQUIREMENTS
        ------------
          antechamber (on PATH)
        """

        # Create temporary working directory and move there.
        old_directory = os.getcwd()
        working_directory = tempfile.mkdtemp()
        os.chdir(working_directory)

        # Write input mol2 file to temporary directory.
        uncharged_molecule_filename = tempfile.mktemp(
            suffix='.mol2',
            dir=working_directory)
        if verbose:
            print "Writing uncharged molecule to %(uncharged_molecule_filename)s" % vars()

        lt.writeMolecule(molecule, uncharged_molecule_filename)

        # Create filename for output mol2 file.
        handle, charged_molecule_filename = tempfile.mkstemp(
            suffix='.mol2', dir=working_directory)

        # Determine net charge of ligand from formal charges.
        formal_charge = lt.formalCharge(molecule)

        # Run antechamber to assign GAFF atom types and charge ligand.
        if netcharge:
            chargestr = '-nc %d' % netcharge
        else:
            chargestr = ''

        command = 'antechamber -i %(uncharged_molecule_filename)s -fi mol2 -o %(charged_molecule_filename)s -fo mol2 -c %(charge_model)s -nc %(formal_charge)d -at gaff %(chargestr)s' % vars(
        )
        if judgetypes:
            command += ' -j %(judgetypes)d' % vars()
        if verbose:
            print command
        output, error = subprocess.Popen(
            command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if verbose:
            print output
        if error:
            raise Exception(error)
        # Read new mol2 file.
        if verbose:
            print "Reading charged molecule from %(charged_molecule_filename)s" % vars()
        charged_molecule = lt.readMolecule(charged_molecule_filename)

        # Clean up temporary working directory.
        if cleanup:
            rm = 'rm -r %s' % working_directory
            subprocess.Popen(
                rm.split(' '),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE).communicate()
        else:
            if verbose:
                print "Work done in %s..." % working_directory

        # Restore old working directory.
        os.chdir(old_directory)

        # Return the charged molecule
        return charged_molecule

    def parametrize(self, **kwargs):
        """
        Parametrize the molecule for amber
        """
        self.prmtop = "molecules/%s.prmtop" % self.name
        self.inpcrd = "molecules/%s.inpcrd" % self.name
        lt.parameterizeForAmber(self.mol, self.prmtop, self.inpcrd, **kwargs)
        return

    def write(self, name=None,wdir="molecules", **kwargs):
        """
        Output the current state of the molecule as a mol2 file.
        """
        if not name:
            name = self.name

        try:
            os.mkdir(wdir)
        except OSError:
            pass

        self.mol2f = "%s/%s.mol2" % (wdir,name)
        lt.writeMolecule(self.mol, self.mol2f, **kwargs)
        return

    def __repr__(self):
        """
        repr function
        """
        return str(vars(self))

    def __str__(self):
        """
        String representation of molecule object
        """

        def _simple(key):
            """
            See if key exists, and show key and value if so.
            """
            try:
                line = "%s : %s\n" % (key, vars(self)[key])
            except:
                line = str()
            return line

        def _fileloc(key):
            """
            See if key exists, resolve path and show abspath to file.
            """
            try:
                line = "%s location: %s\n" % (
                    key, os.path.abspath(vars(self)[key]))
            except:
                line = str()
            return line

        report = str("Reporting molecule \n")
        report_syntax = {'name': _simple('name'),
                         'iupac': _simple('iupac'),
                         'formal_charge': _simple('formal_charge'),
                         'charge_assigned': _simple('charge_assigned'),
                         'mol2f': _fileloc('mol2f'),
                         'prmtop': _fileloc('prmtop'),
                         'inpcrd': _fileloc('inpcrd')
                         }

        for k in vars(self).iterkeys(): #TODO add fixed order
            try:
                report = report + report_syntax[k]
            except KeyError:
                pass

        return report


class MoleculePipeline(object):

    """
    Object for holding molecule data and framework for performing operations
    ARGUMENTS
    ---------
    name - str
        name of the compound
    iupac - str
        IUPAC name of the compound
    formal_charge - int
        the formal charge of the compound
    wdir - str
        working directory
    """

    def __init__(self, name="molecule", iupac="", formal_charge=0, wdir="./"):
        """
        Initialize molecular data structure

        ARGUMENTS
        ---------
            name - str
                name of the compound
            iupac - str
                IUPAC name of the compound
            formal_charge - int
                the formal charge of the compound
            wdir - str
                working directory
        """
        self.name = name
        self.iupac = iupac
        self.formal_charge = int(formal_charge)
        self.charge_assigned = False
        self.wdir = wdir
        self.molstr = list()  # list of molecules as pdb strings

    def add_from_iupac(self, **kwargs):
        """
        Creates a molecule from IUPAC name.

        Requires valid OEChem license

        ARGUMENTS
            **kwargs - keyword args passed directly to createMoleculeFromIUPAC.
        """
        oemol = lt.createMoleculeFromIUPAC(
            self.iupac,
            charge=self.formal_charge,
            **kwargs
            )

        path = self.oemol_write(oemol)
        with open(path, "r") as somefile:
            self.molstr.append(somefile.read())

    def add_with_rdkit(self, filename, filetype, strict=False, **kwargs):
        """Add molecule using rdkit"""
        with open(filename, "r") as inputfile:
            text = inputfile.read()

            if filetype == "auto":
                filetype = os.path.splitext(filename)[1]

            if filetype == "inchi":
                rdmol = Chem.MolFromInchi(text)
            elif filetype == "mol2":
                rdmol = Chem.MolFromMol2File(filename)
            elif filetype == "mol":
                rdmol = Chem.MolFromMolFile(filename)
            elif filetype == "pdb":
                rdmol = Chem.MolFromPdbFile(filename)
            elif filetype in ["smi", "smiles"]:
                rdmol = Chem.MolFromSmiles(text)
            elif filetype == "tpl":
                rdmol = Chem.MolFromTPLFile(filename)
            elif filetype == "smarts":
                if strict:
                    raise IOError("Smarts is pattern, smiles for molecules.")
                else:
                    print "WARNING: Use smiles, ignoring smarts." % filetype
                    return
            else:
                if strict:
                    raise IOError("Filetype (%s) not in rdkit." % filetype)
                else:
                    print "WARNING: Could not filetype %s." % filetype
                    return

        rdmol = Chem.addHs(rdmol)
        self.molstr.append(Chem.MolToPDBBlock)


#  ARGUMENT PARSER

parser = argparse.ArgumentParser(
    description="Transform IUPAC names into molecules with predicated protonation and tautomeric states and get topologies for use with YANK")
inputf = parser.add_mutually_exclusive_group(required=True)

inputf.add_argument(
    "-i",
    dest="file",
    help="Tab separated file with name, iupac name and formal charges of molecules, one molecule per line.",
    metavar="IUPAC",
    type=lambda i: try_opening(
        parser,
         i))
inputf.add_argument(
    "-r",
    dest="filename",
    help="Any molecular input file supported by OEchem.",
    type=str,
    metavar="FILENAME")
inputf.add_argument(
    "-g",
    dest="glob",
    help="Pattern for molecular input file(s) supported by OEchem.",
    type=str,
    metavar="Pattern")

parser.add_argument(
    "-p",
    dest="procedure",
    type=str,
    choices=["sampl4cb7"],
    metavar="proc",
    required=True,
    help="Procedure to put the input through (sampl4cb7 is host-guest).")


"""
MAIN CODE

Plan:

Take IUPAC,mol2, smiles,inchl, etc
|
|
|
|--->Through openeye, convert this to mol2
    |
    |-->Parametrize these molecules through antechamber/parmchk(2) (gaff.mol2)
        |
        |
        |-->frcmod + gaff.mol2 into LEAP -- combine host and guest
            |
            |
            |--> Input prmtop and incrd to Yank

Currently:

Input
    text file with IUPAC names separated by newlines
    or
    molecular filetype, containing 1 molecule each, supported by openeye

Requires
    OpenEye tk
    mmtools.moltools
    Antechamber

To do

    ~ Guest-host integration
    ~ Try using sampl4cb7 set
    ~ Investigate antechamber parametrization bug mentioned by DLM (see lt.assignPartialChargesWithAntechamber)
    ~ Input formats
    ~ Can't split lines by characters used in IUPAC names. Must clearly have TAB separated files.'
    ~ Add more formats
    ~ See how to integrate some parts into Yank in the future

Issues
    Iupac names aren't very robust.
"""


def sampl4cb7(args):
    """The pipeline procedure for Guest-Host"""

    # Extract molecular data from the files
    molecules = list()

    if args.file:
        molstr = args.file.readlines()
        # Cleanup strings.
        molstr = map(str.strip, molstr)
        # Split by tab.
        molstr = [line.split("\t") for line in molstr]
        # Check input file consistency
        for line in molstr:
            if not len(line) == 3:
                raise Exception(
                    "Expecting 3 tab-separated fields per line, verify input.")

        # make a molecular dictionary
        for molecule in molstr:
            molecules.append(OpenEyePipeline(*molecule))

        for mol in molecules:
            mol.create_from_iupac()

    elif args.glob:
        files = glob(args.glob)
        if len(files):
            for f in files:
                molecules.append(OpenEyePipeline())  # initialize as empty
                molecules[-1].create_from_file(args.filename)
        else:
            raise Exception("No files were found.")
    elif args.filename:
        molecules.append(OpenEyePipeline())  # initialize as empty
        molecules[-1].create_from_file(args.filename)
    else:
        raise Exception("Invalid input.")

    # status report.
    for mol in molecules:
        print mol
        mol.write()

    # ensure antechamber is available
    load_amber()

    for mol in molecules:
        oldwd = os.getcwd()
        wd = "molecules/debug"
        try:
            os.mkdir(wd)
        except OSError:
            print "Directory %s exists" % wd

        os.chdir(wd)
        mol.submit_for_pKa()  # dummy function
        mol.gen_conformers(maxconfs=12)
        mol.write_multi("debug")
        mol.load_conformers(pattern="*.mol2")
        print mol.mol
        mol.charge()

#        mol.parametrize(cleanup=False)
#        mol.write()
        os.chdir(oldwd)
        return vars()

if __name__ == "__main__":

    # Process command line input.
    args = parser.parse_args()

    if args.procedure == "sampl4cb7":
        thewholething=sampl4cb7(args)
