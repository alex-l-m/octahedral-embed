Template-based geometry generation for iridium-based phosphorescent emitter molecules, with explicit fac/mer control.
Provides the function `octahedral_embed` to add coordinates to RDKit molecules, using RDKit's distance geometry but constrained to correct octahedral coordination, with control of the isomer.

## Functionality

`octahedral_embed` is intended to be used in place of RDKit's `EmbedMolecule`, solving two problems: that `EmbedMolecule` is not parametrized for metals and will often produce incorrect coordination geometry, and that when it does provide the right coordination, the isomer will be random rather than controlled.
`octahedral_embed` works by matching the input molecule to a template "core", containing the iridium and surrounding atoms.
`octahedral_embed` provides control over the isomer through the `isomer` argument, by selecting a template with the right isomer.

Currently, `octahedral_embed` works for N^C ligands, generating the fac or mer isomers.

Templates are currently based on crystal structures from the Cambridge Structural Database.

## Installation

To install, run in the repo root:

```
pip install .
```

The only dependency is RDKit itself.
This dependency is not enforced, so install into an environment that already has RDKit.

## Usage

### Function `octahedral_embed`

The usage of `octahedral_embed` is similar to RDKit's `EmbedMolecule`, but with an extra `isomer` argument.

```python
octahedral_embed(rdkit_mol, isomer, clearConfs=False)
```

* `rdkit_mol`: RDKit `Mol` (expected to have explicit hydrogens)

* `isomer`: `"fac"` or `"mer"`

* clearConfs: Boolean

Modifies `rdkit_mol` in place, and returns nothing.
Removes all conformers if clearConfs is True, then adds one new conformer.

### Function `ligate`

For convenience, there is also a `ligate` function that will attach ligands to iridium, without adding a geometry.
Ligands must be given as RDKit molecules, with bonds to dummy atoms (symbol "*") where the ligand chelates the metal.

```python
ligate(ligands, metal_atom_element="Ir", metal_atom=None)
```

* `ligands`: list of RDKit `Mol` objects.

* `metal_atom_element`: String with the element of the metal to chelate to.

* `metal_atom`: As an optional alternative to specifying `metal_atom_element`, an RDKit `Atom` object to chelate to.

Returns: an RDKit `Mol` containing the metal with ligands attached, bonding to the metal at the same sites and with the same bond types as they were bonded to the dummy atoms.

## Example

As an example, consider generating the geometry of Ir(ppy)₃.

```python
from rdkit import Chem
from octahedral_embed import ligate, octahedral_embed

# Phenylpyridine ligand
# Use dummy atoms (*) at the chelation sites
# Use a dative bond (->) for the chelating N
ppy = Chem.MolFromSmiles('c1ccc(c(*)c1)c2cccc(n2->*)')

# Ir(ppy)3, with explicit hydrogens to prepare for embedding
irppy3 = Chem.AddHs(ligate([ppy] * 3))

# Save the fac isomer
octahedral_embed(irppy3, isomer='fac')
Chem.MolToMolFile(irppy3, 'irppy3_fac.mol')

# Save the mer isomer
octahedral_embed(irppy3, isomer='mer')
Chem.MolToMolFile(irppy3, 'irppy3_mer.mol')
```
