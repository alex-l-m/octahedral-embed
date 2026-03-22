from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem.rdchem import Mol


def non_iridium_atoms_have_uff_params(mol: Mol) -> bool:
    work = Chem.RWMol(mol)
    iridium_atom_indices = [
        atom.GetIdx()
        for atom in work.GetAtoms()
        if atom.GetSymbol() == "Ir"
    ]
    for atom_idx in reversed(iridium_atom_indices):
        work.RemoveAtom(atom_idx)

    stripped = work.GetMol()
    stripped.UpdatePropertyCache(strict=False)

    if stripped.GetNumAtoms() == 0:
        return True

    return rdForceFieldHelpers.UFFHasAllMoleculeParams(stripped)
