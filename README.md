Template-based geometry generation for iridium-based phosphorescent emitter molecules, with explicit fac/mer control.
Provides the function `octahedral_embed` to add coordinates to RDKit molecules, using RDKit's distance geometry but constrained to correct octahedral coordination, with control of the isomer.

## Functionality

`octahedral_embed` is intended to be used in place of RDKit's `EmbedMolecule`.
It works by matching the input molecule to a template "core", containing the iridium and surrounding atoms.
It provides control over the isomer by selecting an isomer-specific template.
The templates were extracted from crystal structures from the Cambridge Structural Database.

Currently, `octahedral_embed` can generate tris-bidentate complexes with three N^C ligands or three carbene ligands.
The `isomer` argument controls whether the output is the fac or mer isomer.

`octahedral_embed` will become obsolete when (if?) the RDKit properly implements non-tetrahedral stereochemistry and includes metals in force fields.

## Installation

To install, run in the repo root:

```
pip install .
```

The package depends on RDKit.

## Usage

### Function `octahedral_embed`

The usage of `octahedral_embed` is similar to RDKit's `EmbedMolecule`, but with an extra `isomer` argument.

```python
octahedral_embed(rdkit_mol, isomer, clearConfs=True)
```

* `rdkit_mol`: RDKit `Mol` (expected to have explicit hydrogens)

* `isomer`: `"fac"` or `"mer"`

* clearConfs: Boolean

Modifies `rdkit_mol` in place.
Removes all conformers if clearConfs is True, then adds one new conformer.

Returns the index of the added conformer.

### Molecule-construction utilities

Utilities for building input molecules live in `octahedral_embed.construction`.
These currently include `ligate` and `make_bonds_dative`.

#### Function `ligate`

For convenience, there is also a `ligate` function that will attach ligands to iridium, without adding a geometry.
Ligands must be given as RDKit molecules, with bonds to dummy atoms (symbol `"*"`) where the ligand chelates the metal.

```python
ligate(ligands, metal_atom_element="Ir", metal_atom=None)
```

* `ligands`: list of RDKit `Mol` objects.

* `metal_atom_element`: String with the element of the metal to chelate to.

* `metal_atom`: As an optional alternative to specifying `metal_atom_element`, an RDKit `Mol` containing the metal atom to chelate to.

Returns: an RDKit `Mol` containing the metal with ligands attached, bonding to the metal at the same sites and with the same bond types as they were bonded to the dummy atoms.

## Example

As an example, consider generating the geometry of Ir(ppy)₃.

```python
from rdkit import Chem
from octahedral_embed import octahedral_embed
from octahedral_embed.construction import ligate

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

For a carbene emitter, here is a simple example for Ir(pmb)₃.

```python
from rdkit import Chem
from octahedral_embed import octahedral_embed
from octahedral_embed.construction import ligate

# pmb ligand
# Use dummy atoms (*) at the chelation sites
# Use a double bond (=*) for the carbene carbon
pmb = Chem.MolFromSmiles('*c1ccccc1N1C(=*)N(C)c2ccccc21')

irpmb3 = Chem.AddHs(ligate([Chem.Mol(pmb) for _ in range(3)]))

octahedral_embed(irpmb3, isomer='fac')
Chem.MolToMolFile(irpmb3, 'irpmb3_fac.mol')
```
