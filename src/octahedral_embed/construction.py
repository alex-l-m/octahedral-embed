from functools import reduce
from typing import Optional, Sequence

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Mol, RWMol
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import CombineMols, Kekulize, SanitizeMol

def set_single_iridium_formal_charge_zero(mol: Chem.Mol) -> None:
    (ir_idx,) = (atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "Ir")
    mol.GetAtomWithIdx(ir_idx).SetFormalCharge(0)
    mol.UpdatePropertyCache(strict=False)


def ligate(
    ligands: Sequence[Mol],
    metal_atom_element: str = "Ir",
    metal_atom: Optional[Mol] = None,
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
        if (
            bond.GetBeginAtom().GetSymbol() == target_elem
            and bond.GetEndAtom().GetSymbol() in ["N", "P"]
            and bond.GetEndAtom().GetFormalCharge() == 1
        ):
            iridium = bond.GetBeginAtom()
            nitrogen = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif (
            bond.GetEndAtom().GetSymbol() == target_elem
            and bond.GetBeginAtom().GetSymbol() in ["N", "P"]
            and bond.GetBeginAtom().GetFormalCharge() == 1
        ):
            iridium = bond.GetEndAtom()
            nitrogen = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
        if (
            bond.GetBeginAtom().GetSymbol() == target_elem
            and bond.GetEndAtom().GetSymbol() == "C"
            and bond.GetEndAtom().GetTotalValence() == 3
        ):
            iridium = bond.GetBeginAtom()
            carbene = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif (
            bond.GetEndAtom().GetSymbol() == target_elem
            and bond.GetBeginAtom().GetSymbol() == "C"
            and bond.GetBeginAtom().GetTotalValence() == 3
        ):
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
