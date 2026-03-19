from collections.abc import Iterable

from rdkit.Chem.rdchem import BondStereo, ChiralType, Mol

_ATOM_STEREO_PROPS = (
    "_CIPCode",
    "_CIPRank",
    "_ChiralityPossible",
    "_chiralPermutation",
)
_BOND_STEREO_PROPS = ("_CIPCode",)


def _clear_atom_stereo(atom) -> None:
    atom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
    for prop in _ATOM_STEREO_PROPS:
        if atom.HasProp(prop):
            atom.ClearProp(prop)


def _clear_bond_stereo(bond) -> None:
    bond.SetStereo(BondStereo.STEREONONE)
    for prop in _BOND_STEREO_PROPS:
        if bond.HasProp(prop):
            bond.ClearProp(prop)


def clear_stereo_on_atom_indices(mol: Mol, atom_indices: Iterable[int]) -> None:
    atom_index_set = set(atom_indices)
    for atom_idx in atom_index_set:
        _clear_atom_stereo(mol.GetAtomWithIdx(atom_idx))
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtomIdx() in atom_index_set
            and bond.GetEndAtomIdx() in atom_index_set
        ):
            _clear_bond_stereo(bond)


def clear_stereo_on_element(mol: Mol, element: str) -> None:
    atom_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == element
    ]
    clear_stereo_on_atom_indices(mol, atom_indices)
