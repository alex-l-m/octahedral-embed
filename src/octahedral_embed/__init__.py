import os
from typing import Any
from functools import reduce
from typing import List, Optional, Sequence, Literal, cast
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, Conformer, Mol, Atom
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolFromSmarts, MolFromSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.rdmolops import RemoveStereochemistry, CombineMols, SanitizeMol, Kekulize
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from rdkit.Chem.rdMolAlign import AlignMol

# NOTE:
# In some RDKit builds, UFFGetMoleculeForceField can return a
# ForceFields::PyForceField instance. The Boost.Python class registration for
# that type lives in rdkit.ForceField.rdForceField, and without importing it
# first, you can see:
#   TypeError: No Python class registered for C++ class ForceFields::PyForceField
# Importing rdkit.ForceField.rdForceField ensures the wrapper classes are
# registered before we call into force-field code.
import rdkit.ForceField.rdForceField  # noqa: F401


def ConstrainedEmbed_withParams(
    mol: Mol,
    core: Mol,
    useTethers: bool = True,
    coreConfId: int = -1,
    randomseed: Optional[int] = None,
    getForceField: Any = rdForceFieldHelpers.UFFGetMoleculeForceField,
    params: Optional[rdDistGeom.EmbedParameters] = None,
) -> Mol:
    """
    Constrained embedding that honors an EmbedParameters object (timeout, maxIterations, ...).

    Parameters
    ----------
    mol, core : RDKit Mol
        'core' must be a substructure of 'mol'.
    params : rdDistGeom.EmbedParameters, optional
        e.g. rdDistGeom.ETKDGv3(). If None, defaults to rdDistGeom.ETKDGv3().

    Returns
    -------
    mol : RDKit Mol
        With a single embedded conformer and property 'EmbedRMS'.
    """
    match = mol.GetSubstructMatch(core)
    if not match:
        raise ValueError("molecule doesn't match the core")

    # Build the coordMap using the core conformation we were asked to use:
    coreConf = core.GetConformer(coreConfId)
    coordMap = {mol_idx: coreConf.GetAtomPosition(core_atom_idx)
                for core_atom_idx, mol_idx in enumerate(match)}

    # Start from user params (or a sensible default)
    if params is None:
        params = rdDistGeom.ETKDGv3()
    p = params

    # Seed + constraints go onto the cloned params:
    if randomseed is not None:
        p.randomSeed = int(randomseed)
    p.SetCoordMap(coordMap)

    # Embed using the parameter object overload:
    confId = rdDistGeom.EmbedMolecule(mol, p)
    if confId < 0:
        raise ValueError("Could not embed molecule.")

    # Map for alignment: (probeAtomId, refAtomId)
    algMap = [(mol_idx, core_atom_idx) for core_atom_idx, mol_idx in enumerate(match)]

    if not useTethers:
        # Distance-constraint cleanup on the newly generated conformer:
        ff = getForceField(mol, confId=confId)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.0)

        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1

        rms = AlignMol(mol, core, atomMap=algMap, prbCid=confId, refCid=coreConfId)
    else:
        # Align first, then add "tethers" to the core atom positions:
        rms = AlignMol(mol, core, atomMap=algMap, prbCid=confId, refCid=coreConfId)

        ff = getForceField(mol, confId=confId)
        for core_atom_idx in range(core.GetNumAtoms()):
            pt = coreConf.GetAtomPosition(core_atom_idx)
            pIdx = ff.AddExtraPoint(pt.x, pt.y, pt.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[core_atom_idx], 0, 0, 100.0)

        ff.Initialize()
        n = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1

        rms = AlignMol(mol, core, atomMap=algMap, prbCid=confId, refCid=coreConfId)

    mol.SetProp("EmbedRMS", str(rms))
    return mol

def ligate(
    ligands: Sequence[Mol],
    metal_atom_element: str = "Ir",
    metal_atom : Optional[Mol] = None,
) -> Mol:
    for ligand in ligands:
        ligand.RemoveAllConformers()
        Kekulize(ligand)
    if metal_atom is None:
        metal_atom = MolFromSmiles(f"[{metal_atom_element}]")
    assert metal_atom is not None
    # Create a molecule that contains the metal atom as well as all the
    # ligands, but without bonds between the metal and the ligands
    mol = reduce(CombineMols, ligands, metal_atom)
    # Record the index of the metal atom, so that bonds can be created to it
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    # This will not work if the ligands also contain the metal atom!
    metal_atom_index = elements.index(metal_atom_element)
    # Create an editable molecule and begin batch editing, which makes
    # functions available for adding and removing bonds as well as removing the
    # dummy atoms, and should keep the indexes stable as we do so
    editable = Chem.EditableMol(mol)
    editable.BeginBatchEdit()
    # Check each bond to see if it is a coordination site, and if so, bond the
    # atom to the metal
    for bond in mol.GetBonds():
        # Figure out which, if any, of the atoms in the bond is the coordinating atom
        if bond.GetBeginAtom().GetSymbol() == "*":
            dummy_atom_index = bond.GetBeginAtomIdx()
            coordination_atom_index = bond.GetEndAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == "*":
            dummy_atom_index = bond.GetEndAtomIdx()
            coordination_atom_index = bond.GetBeginAtomIdx()
        else:
            continue
        # Add a bond to metal of the same type as the bond to the dummy atom
        # Has to have metal second for the dative bonds to work
        # Order matters for dative bonds
        # https://www.rdkit.org/docs/RDKit_Book.html#dative-bonds
        editable.AddBond(coordination_atom_index, metal_atom_index, bond.GetBondType())
        # Remove the bond between the coordinating atom and the dummy atom
        editable.RemoveBond(dummy_atom_index, coordination_atom_index)
    # Remove all the dummy atoms
    for i, element in enumerate(elements):
        if element == "*":
            editable.RemoveAtom(i)
    # Apply the changes and get the molecule
    editable.CommitBatchEdit()
    outmol = editable.GetMol()
    # Shouldn't have any conformers but just in case
    outmol.RemoveAllConformers()
    # Probably already done by GetMol but just in case
    SanitizeMol(outmol)
    return outmol

def make_bonds_dative(mol: Mol, target_elem: str = "Ir") -> Mol:
    editable_mol = RWMol(mol)

    # If you don't make a list, it loops infinitely over the bonds it's creating
    for bond in list(editable_mol.GetBonds()):
        iridium: Optional[Atom] = None
        nitrogen: Optional[Atom] = None
        carbene: Optional[Atom] = None
        start_idx: Optional[int] = None
        end_idx: Optional[int] = None
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() in ["N", "P"] and \
                bond.GetEndAtom().GetFormalCharge() == 1:
            iridium = bond.GetBeginAtom()
            nitrogen = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() in ["N", "P"] and \
                bond.GetBeginAtom().GetFormalCharge() == 1:
            iridium = bond.GetEndAtom()
            nitrogen = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() == "C" and \
                bond.GetEndAtom().GetTotalValence() == 3:
            iridium = bond.GetBeginAtom()
            carbene = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() == "C" and \
                bond.GetBeginAtom().GetTotalValence() == 3:
            iridium = bond.GetEndAtom()
            carbene = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

        # Convert from CSD bond conventions to RDKit-friendly ones.
        #
        # * N+/P+ - Ir single bonds are best represented as dative (N/P -> Ir)
        #   to respect RDKit valence rules.
        # * Carbene C - Ir single bonds (carbon total valence == 3) are encoded
        #   here as *double bonds* (C=Ir) rather than dative bonds.
        if nitrogen is not None:
            # Replace N+ - Ir with N -> Ir
            nitrogen.SetFormalCharge(0)
            assert start_idx is not None and end_idx is not None
            editable_mol.RemoveBond(start_idx, end_idx)
            editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DATIVE)
        elif iridium is not None and carbene is not None:
            assert start_idx is not None and end_idx is not None
            editable_mol.RemoveBond(start_idx, end_idx)
            editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DOUBLE)

    outmol = editable_mol.GetMol()
    Chem.SanitizeMol(outmol)

    return outmol

def load_template(filename: str) -> Mol:
    inpath = os.path.join(__path__[0], filename)
    # Load, sanitized, without hydrogens
    # Loading without hydrogens is fine because they're not going to be used
    # for matching the template
    # Sanitizing is fine, even though these templates are from the CSD, because
    # while not all CSD molecules sanitize, I'm only using the ones that do as
    # templates
    raw_mol = MolFromMol2File(inpath)
    assert raw_mol is not None
    # Stereochemistry makes template matching more strict
    RemoveStereochemistry(raw_mol)
    # Convert from CSD bond conventions, to respect RDKit valence rules
    dative_mol = make_bonds_dative(raw_mol)
    # Center the molecule
    # First, retrieve the index of the iridium atom
    target_index: Optional[int] = None
    for atom in dative_mol.GetAtoms():
        element = atom.GetSymbol()
        if element == 'Ir':
            target_index = atom.GetIdx()
            break
    assert target_index is not None
    # Then, retrieve the coordinates of the iridium atom in the geometry
    geometry = dative_mol.GetConformer()
    target_point = geometry.GetAtomPosition(target_index)
    # Finally, use these coordinates as the center
    CanonicalizeConformer(geometry, center=target_point)
    return dative_mol

def transfer_conformation(mol: Mol, substruct: Mol, conformer: int = 0) -> None:
    '''Given a molecule, and a second molecule which is a substructure of the
    first, assign coordinates to the substructure based on the matching part of
    the original molecule'''
    match = mol.GetSubstructMatch(substruct)
    substruct_conformation = Conformer(substruct.GetNumAtoms())
    for i, index in enumerate(match):
        point = mol.GetConformer(conformer).GetAtomPosition(index)
        substruct_conformation.SetAtomPosition(i, point)
    substruct.AddConformer(substruct_conformation)

fac = load_template("OHUZEW.mol2")
mer = load_template("OHUZIA.mol2")

carbene_fac = load_template("MAXYIU.mol2")
carbene_mer = load_template("MAXYOA.mol2")

facmer_skeleton_smarts = (
    "[Ir]123"
    "(~[n]~[a]~[a]~[c]~1)"
    "(~[n]~[a]~[a]~[c]~2)"
    "(~[n]~[a]~[a]~[c]~3)"
)

fac_skeleton = cast(Mol, MolFromSmarts(facmer_skeleton_smarts))
transfer_conformation(fac, fac_skeleton)

mer_skeleton = cast(Mol, MolFromSmarts(facmer_skeleton_smarts))
transfer_conformation(mer, mer_skeleton)

carbene_skeleton_smarts = "[Ir]135(=[CH0](~N(~*)~*~2)~N(~*~2)~c~c~1)(=[CH0](~N(~*)~*~4)~N(~*~4)~c~c~3)(=[CH0](~N(~*)~*~6)~N(~*~6)~c~c~5)"
carbene_fac_skeleton = cast(Mol, MolFromSmarts(carbene_skeleton_smarts))
transfer_conformation(carbene_fac, carbene_fac_skeleton)
carbene_mer_skeleton = cast(Mol, MolFromSmarts(carbene_skeleton_smarts))
transfer_conformation(carbene_mer, carbene_mer_skeleton)

fac_skeletons: List[Mol] = [fac_skeleton, carbene_fac_skeleton]
mer_skeletons: List[Mol] = [mer_skeleton, carbene_mer_skeleton]

def octahedral_embed(
    mol: Mol,
    isomer: Literal['fac', 'mer'],
    clearConfs: bool = True,
) -> int:
    # Make a copy of the molecule to avoid side effects (in particular, removing stereochemistry)
    work = Mol(mol)

    # Needed for some of the mol2 files I got from CSD
    # Will not be able to embed with stereochemistry
    RemoveStereochemistry(work)

    work.RemoveAllConformers()

    if isomer == "fac":
        skeletons = fac_skeletons
    elif isomer == "mer":
        skeletons = mer_skeletons
    else:
        raise ValueError(f"Isomer should be \"mer\" or \"fac\", given {isomer}")
    finished = False
    for skeleton in skeletons:
        if len(work.GetSubstructMatch(skeleton)) > 0:
            ps = rdDistGeom.ETKDGv3()

            # Carbene embedding with a large template gives output "Could not
            # triangle bounds smooth molecule" and raises a ValueError. But
            # with a small template the imidazole is horribly twisted, probably
            # because it thinks the atoms are aliphatic. Ignoring smoothing
            # failures with the large template, it works. Now I'm hoping that
            # if I encode carbenes more reasonably, I won't have these issues.
            # However, ignoring smoothing failures is always an option later if
            # required.
            #ps.ignoreSmoothingFailures = True
            
            ps.timeout = 5
            ps.maxIterations = 20

            work = ConstrainedEmbed_withParams(work, skeleton, params=ps)
            finished = True
            break
    if not finished:
        raise ValueError("Doesn't match templates")
    # Copy the conformer back to the original molecule
    new_conf = work.GetConformer()
    if clearConfs:
        mol.RemoveAllConformers()
    return mol.AddConformer(new_conf, assignId=True)
