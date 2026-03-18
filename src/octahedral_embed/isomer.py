from itertools import combinations
from typing import Literal, Optional, Sequence, Tuple

from rdkit.Chem.rdMolTransforms import GetAngleDeg
from rdkit.Chem.rdchem import Mol

from .templates import FAMILY_SPECS, FamilySpec

# Simple fac/mer discriminator for equivalent donor atoms in an octahedral
# tris-bidentate complex. Ideal fac gives three ~90° same-donor angles;
# ideal mer gives one ~180° angle. 135° is the midpoint and works well for
# distorted structures.
TRANS_ANGLE_CUTOFF_DEG = 135.0


def find_supported_family(mol: Mol) -> Tuple[Optional[FamilySpec], Tuple[int, ...]]:
    for spec in FAMILY_SPECS:
        match = mol.GetSubstructMatch(spec.query)
        if match:
            return spec, tuple(match)
    return None, tuple()



def pairwise_angles_at_center(
    mol: Mol,
    center_idx: int,
    atom_indices: Sequence[int],
    conf_id: int = 0,
) -> list[float]:
    conf = mol.GetConformer(conf_id)
    return sorted(
        float(GetAngleDeg(conf, i, center_idx, j))
        for i, j in combinations(atom_indices, 2)
    )



def fac_or_mer(mol: Mol, conf_id: int = 0) -> Literal["fac", "mer"]:
    spec, match = find_supported_family(mol)
    if spec is None:
        raise ValueError("Molecule does not match a supported fac/mer family")

    basis_indices = spec.donor_sets[spec.fac_mer_basis]
    basis_atom_indices = [match[i] for i in basis_indices]
    basis_angles = pairwise_angles_at_center(
        mol,
        center_idx=match[0],
        atom_indices=basis_atom_indices,
        conf_id=conf_id,
    )
    return "mer" if max(basis_angles) >= TRANS_ANGLE_CUTOFF_DEG else "fac"



def assert_isomer(
    mol: Mol,
    expected_isomer: Literal["fac", "mer"],
    conf_id: int = 0,
) -> None:
    observed_isomer = fac_or_mer(mol, conf_id=conf_id)
    if observed_isomer != expected_isomer:
        raise ValueError(
            f'Embedded molecule has {observed_isomer} geometry, expected {expected_isomer}'
        )
