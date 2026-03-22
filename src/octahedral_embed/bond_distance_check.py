from rdkit.Chem import rdForceFieldHelpers, rdMolTransforms
from rdkit.Chem.rdchem import Mol

DEFAULT_MAX_UFF_RELATIVE_BOND_STRETCH = 0.30


def bond_distance_check(
    mol: Mol,
    confId: int = -1,
    max_relative_stretch: float = DEFAULT_MAX_UFF_RELATIVE_BOND_STRETCH,
) -> bool:
    """Return True when all checked bonds are close to their UFF rest lengths.

    For each non-hydrogen bond with available UFF stretch parameters, compare
    the bond length to the corresponding UFF r0. Bonds are accepted when |d -
    r0| / r0 <= max_relative_stretch.
    """

    conf = mol.GetConformer(confId)

    for bond in mol.GetBonds():

        params = rdForceFieldHelpers.GetUFFBondStretchParams(
            mol,
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
        )
        if params is None:
            continue

        _, r0 = params
        distance = rdMolTransforms.GetBondLength(
            conf,
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
        )
        if abs(distance - r0) / r0 > max_relative_stretch:
            return False

    return True
