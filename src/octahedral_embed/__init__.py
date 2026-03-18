from typing import Literal

from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import RemoveStereochemistry

from .constrained_embed import ConstrainedEmbed_withParams
from .isomer import assert_isomer
from .templates import (
    fac_skeletons,
    mer_skeletons,
)
from .construction import set_single_iridium_formal_charge_zero

def octahedral_embed(
    mol: Mol,
    isomer: Literal['fac', 'mer'],
    clearConfs: bool = True,
) -> int:
    # Make a copy of the molecule to avoid side effects (in particular, removing stereochemistry)
    work = Mol(mol)

    # Set iridium formal charge to zero
    # RDKit assigns some crazy formal charges to iridium
    # This decreases the rate of failures when testing on molecules
    # loaded from the CSD
    set_single_iridium_formal_charge_zero(work)

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
            assert_isomer(work, isomer)
            finished = True
            break
    if not finished:
        raise ValueError("Doesn't match templates")
    # Copy the conformer back to the original molecule
    new_conf = work.GetConformer()
    if clearConfs:
        mol.RemoveAllConformers()
    return mol.AddConformer(new_conf, assignId=True)
