from typing import Any, Optional

from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem.rdchem import Mol, SubstructMatchParameters

from .isoctahedral import isoctahedral

# NOTE:
# In some RDKit builds, UFFGetMoleculeForceField can return a
# ForceFields::PyForceField instance. The Boost.Python class registration for
# that type lives in rdkit.ForceField.rdForceField, and without importing it
# first, you can see:
#   TypeError: No Python class registered for C++ class ForceFields::PyForceField
# Importing rdkit.ForceField.rdForceField ensures the wrapper classes are
# registered before we call into force-field code.
import rdkit.ForceField.rdForceField


def ConstrainedEmbed_withParams(
    mol: Mol,
    core: Mol,
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
    # Don't use chirality in the substruct match
    # Hopefully my code stripping stereochemistry from the iridium atom is no
    # longer necessary to achieve a match (though it may still be necessary for
    # embedding)
    match_params = SubstructMatchParameters()
    match_params.useChirality = False

    match = mol.GetSubstructMatch(core, match_params)
    if not match:
        raise ValueError("molecule doesn't match the core")

    # Build the coordMap using the core conformation we were asked to use:
    coreConf = core.GetConformer(coreConfId)
    coordMap = {
        mol_idx: coreConf.GetAtomPosition(core_atom_idx)
        for core_atom_idx, mol_idx in enumerate(match)
    }

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

    if more != 0:
        raise ValueError('Could not optimize molecule')

    if not isoctahedral(mol, confId):
        raise ValueError('Embedded molecule does not satisfy octahedral geometry')

    rms = AlignMol(mol, core, atomMap=algMap, prbCid=confId, refCid=coreConfId)

    mol.SetProp("EmbedRMS", str(rms))
    return mol
