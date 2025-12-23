Geometry generation for iridium-based phosphorescent emitter molecules.
Provides the function `octahedral_embed` to add coordinates to RDKit molecules, using RDKit's distance geometry, with control of the isomer.

## Functionality

`octahedral_embed` is intended to be used in place of RDKit's `EmbedMolecule`, solving two problems: that `EmbedMolecule` is not parametrized for metals and will often produce incorrect coordination geometry, and that when it does provide the right coordination, the isomer will be random rather than controlled.
`octahedral\_embed` works by matching the input molecule to a template "core", containing the iridium and surrounding atoms, extracted from a crystal structure.
`octahedral\_embed` provides control over the isomer through the `isomer` argument, by selecting a template with the right isomer.

Currently, octahedral\_embed works for NC ligands, generating in the fac or mer isomers.

There are no dependencies besides RDKit itself.

## Usage

### Function `octahedral_embed`

The usage of `octahedral_embed` is similar to RDKit's `EmbedMolecule`, but with an extra `isomer` argument.

```
octahedral_embed(rdkit_mol, isomer)
```

* `mol`: RDKit `Mol` (expected to have explicit hydrogens)

* `isomer`: `"fac"` or `"mer"`

Adds a conformer in-place, replacing any existing conformers.
No return value.

### Function `ligate`

For convenience, there is also a `ligate` function that will attach ligands to iridium, without adding a geometry.
Ligands must be given as RDKit molecules, with bonds to dummy atoms (symbol "*") where the ligand chelates the metal.

```
ligate(ligands, metal_atom_element="Ir", metal_atom=None)
```

* `ligands`: list of RDKit `Mol` objects.

* `metal_atom_element`: String with the element of the metal to chelate to.

* `metal_atom`: As an optional alternative to specifying `metal_atom_element`, an RDKit `Atom` object to chelate to.

Returns: an RDKit `Mol` containing the metal with ligands attached, bonding to the metal at the same sites and with the same bond types as they were bonded to the dummy atoms.

## Example

As an example, consider generating the geometry of Ir(ppy)₃.

```
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
