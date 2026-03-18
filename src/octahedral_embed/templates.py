import os
from dataclasses import dataclass
from typing import Dict, Tuple, cast

from rdkit.Chem.rdchem import Conformer, Mol
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolFromSmarts
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdmolops import RemoveStereochemistry

from .construction import make_bonds_dative, set_single_iridium_formal_charge_zero

FACMER_SKELETON_SMARTS = (
    "[Ir]123"
    "(~[n]~[a]~[a]~[c]~1)"
    "(~[n]~[a]~[a]~[c]~2)"
    "(~[n]~[a]~[a]~[c]~3)"
)

CARBENE_SKELETON_SMARTS = (
    "[Ir]135"
    "(=[CH0](~N(~*)~*~2)~N(~*~2)~c~c~1)"
    "(=[CH0](~N(~*)~*~4)~N(~*~4)~c~c~3)"
    "(=[CH0](~N(~*)~*~6)~N(~*~6)~c~c~5)"
)


@dataclass(frozen=True)
class FamilySpec:
    name: str
    query: Mol
    donor_sets: Dict[str, Tuple[int, int, int]]
    fac_mer_basis: str


def _query_from_smarts(smarts: str) -> Mol:
    query = MolFromSmarts(smarts)
    if query is None:
        raise ValueError(f"Could not parse SMARTS: {smarts}")
    return cast(Mol, query)


FAMILY_SPECS = (
    FamilySpec(
        name="N^C",
        query=_query_from_smarts(FACMER_SKELETON_SMARTS),
        donor_sets={
            "N": (1, 5, 9),
            "cyclometalated_C": (4, 8, 12),
        },
        fac_mer_basis="N",
    ),
    FamilySpec(
        name="carbene",
        query=_query_from_smarts(CARBENE_SKELETON_SMARTS),
        donor_sets={
            "carbene_C": (1, 9, 17),
            "cyclometalated_C": (8, 16, 24),
        },
        fac_mer_basis="carbene_C",
    ),
)


def load_template(filename: str) -> Mol:
    inpath = os.path.join(os.path.dirname(__file__), filename)
    # Load, sanitized, without hydrogens
    # Loading without hydrogens is fine because they're not going to be used
    # for matching the template
    # Sanitizing is fine, even though these templates are from the CSD, because
    # while not all CSD molecules sanitize, I'm only using the ones that do as
    # templates
    raw_mol = MolFromMol2File(inpath)
    assert raw_mol is not None
    # Set iridium formal charge to zero
    set_single_iridium_formal_charge_zero(raw_mol)
    # Stereochemistry makes template matching more strict
    RemoveStereochemistry(raw_mol)
    # Convert from CSD bond conventions, to respect RDKit valence rules
    dative_mol = make_bonds_dative(raw_mol)
    # Center the molecule
    # First, retrieve the index of the iridium atom
    target_index = None
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



def load_skeleton(template: Mol, smarts: str) -> Mol:
    skeleton = cast(Mol, MolFromSmarts(smarts))
    assert skeleton is not None
    transfer_conformation(template, skeleton)
    return skeleton


fac = load_template("OHUZEW.mol2")
mer = load_template("OHUZIA.mol2")

carbene_fac = load_template("MAXYIU.mol2")
carbene_mer = load_template("MAXYOA.mol2")

fac_skeleton = load_skeleton(fac, FACMER_SKELETON_SMARTS)
mer_skeleton = load_skeleton(mer, FACMER_SKELETON_SMARTS)
carbene_fac_skeleton = load_skeleton(carbene_fac, CARBENE_SKELETON_SMARTS)
carbene_mer_skeleton = load_skeleton(carbene_mer, CARBENE_SKELETON_SMARTS)

fac_skeletons = [fac_skeleton, carbene_fac_skeleton]
mer_skeletons = [mer_skeleton, carbene_mer_skeleton]
