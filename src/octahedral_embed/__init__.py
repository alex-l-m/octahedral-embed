from typing import Literal

from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdchem import Mol, SubstructMatchParameters

from .constrained_embed import ConstrainedEmbed_withParams
from .isomer import assert_isomer
from .templates import (
    fac_skeletons,
    mer_skeletons,
)
from .construction import set_single_iridium_formal_charge_zero
from .stereo import clear_stereo_on_element
from .input_checks import non_iridium_atoms_have_uff_params

def octahedral_embed(
    mol: Mol,
    isomer: Literal['fac', 'mer'],
    clearConfs: bool = True,
) -> int:

    # Check input
    if not non_iridium_atoms_have_uff_params(mol):
        raise ValueError("Molecule not supported by UFF")

    # Make a copy of the molecule to avoid side effects (in particular, removing stereochemistry)
    work = Mol(mol)

    # Set iridium formal charge to zero
    # RDKit assigns some crazy formal charges to iridium
    # This decreases the rate of failures when testing on molecules
    # loaded from the CSD
    set_single_iridium_formal_charge_zero(work)

    # Needed for some of the mol2 files I got from CSD
    # Will not be able to embed with stereochemistry
    # Only clearing on iridium seems to get most or all of the benefit
    clear_stereo_on_element(work, "Ir")

    work.RemoveAllConformers()

    if isomer == "fac":
        skeletons = fac_skeletons
    elif isomer == "mer":
        skeletons = mer_skeletons
    else:
        raise ValueError(f"Isomer should be \"mer\" or \"fac\", given {isomer}")
    finished = False

    # Don't use chirality in the substruct match
    # Hopefully my code stripping stereochemistry from the iridium atom is no
    # longer necessary to achieve a match (though it may still be necessary for
    # embedding)
    match_params = SubstructMatchParameters()
    match_params.useChirality = False

    for skeleton in skeletons:
        if len(work.GetSubstructMatch(skeleton, params=match_params)) > 0:
            ps = rdDistGeom.ETKDGv3()

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
