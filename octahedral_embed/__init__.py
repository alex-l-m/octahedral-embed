import os.path
from functools import reduce
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, Conformer, Mol
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolFromSmarts, MolFromSmiles
from rdkit.Chem.AllChem import ConstrainedEmbed
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.rdmolops import RemoveStereochemistry, CombineMols, SanitizeMol, Kekulize
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer

def ligate(ligands, metal_atom_element = "Ir", metal_atom = None):
    for ligand in ligands:
        ligand.RemoveAllConformers()
        Kekulize(ligand)
    if metal_atom is None:
        metal_atom = MolFromSmiles(f"[{metal_atom_element}]")
    # Create a molecule that contains the metal atom as well as all the
    # ligands, but without bonds between the metal and the ligands
    mol = reduce(CombineMols, ligands, metal_atom)
    # Record the index of the metal atom, so that bonds can be created to it
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    # This will not work if the ligands also contain the metal atom!
    metal_atom_index = elements.index(metal_atom_element)
    # Create an editable molecule and begin batch editing, which makes
    # functions available for adding and removing bonds as well as removing the
    # dummy atoms, and should keep the indexes stable as we do so
    editable = Chem.EditableMol(mol)
    editable.BeginBatchEdit()
    # Check each bond to see if it is a coordination site, and if so, bond the
    # atom to the metal
    for bond in mol.GetBonds():
        # Figure out which, if any, of the atoms in the bond is the coordinating atom
        if bond.GetBeginAtom().GetSymbol() == "*":
            dummy_atom_index = bond.GetBeginAtomIdx()
            coordination_atom_index = bond.GetEndAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == "*":
            dummy_atom_index = bond.GetEndAtomIdx()
            coordination_atom_index = bond.GetBeginAtomIdx()
        else:
            continue
        # Add a bond to metal of the same type as the bond to the dummy atom
        # Has to have metal second for the dative bonds to work
        # Order matters for dative bonds
        # https://www.rdkit.org/docs/RDKit_Book.html#dative-bonds
        editable.AddBond(coordination_atom_index, metal_atom_index, bond.GetBondType())
        # Remove the bond between the coordinating atom and the dummy atom
        editable.RemoveBond(dummy_atom_index, coordination_atom_index)
    # Remove all the dummy atoms
    for i, element in enumerate(elements):
        if element == "*":
            editable.RemoveAtom(i)
    # Apply the changes and get the molecule
    editable.CommitBatchEdit()
    outmol = editable.GetMol()
    # Shouldn't have any conformers but just in case
    outmol.RemoveAllConformers()
    # Probably already done by GetMol but just in case
    SanitizeMol(outmol)
    return outmol

def make_bonds_dative(mol, target_elem="Ir"):
    editable_mol = RWMol(mol)

    # If you don't make a list, it loops infinitely over the bonds it's creating
    for bond in list(editable_mol.GetBonds()):
        iridium = None
        nitrogen = None
        carbene = None
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() in ["N", "P"] and \
                bond.GetEndAtom().GetFormalCharge() == 1:
            iridium = bond.GetBeginAtom()
            nitrogen = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() in ["N", "P"] and \
                bond.GetBeginAtom().GetFormalCharge() == 1:
            iridium = bond.GetEndAtom()
            nitrogen = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() == "C" and \
                bond.GetEndAtom().GetTotalValence() == 3:
            iridium = bond.GetBeginAtom()
            carbene = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() == "C" and \
                bond.GetBeginAtom().GetTotalValence() == 3:
            iridium = bond.GetEndAtom()
            carbene = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

        if nitrogen is not None:
            # Replace N+ - Ir with N -> Ir
            nitrogen.SetFormalCharge(0)

        if iridium is not None and (nitrogen is not None or carbene is not None):
            editable_mol.RemoveBond(start_idx, end_idx)
            editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DATIVE)

    outmol = editable_mol.GetMol()
    Chem.SanitizeMol(outmol)

    return outmol

def load_template(filename):
    inpath = os.path.join(__path__[0], filename)
    # Load, sanitized, without hydrogens
    # Loading without hydrogens is fine because they're not going to be used
    # for matching the template
    # Sanitizing is fine, even though these templates are from the CSD, because
    # while not all CSD molecules sanitize, I'm only using the ones that do as
    # templates
    raw_mol = MolFromMol2File(inpath)
    # Stereochemistry makes template matching more strict
    RemoveStereochemistry(raw_mol)
    # Convert from CSD bond conventions, to respect RDKit valence rules
    dative_mol = make_bonds_dative(raw_mol)
    # Center the molecule
    # First, retrieve the index of the iridium atom
    target_index = None
    for atom in dative_mol.GetAtoms():
        element = atom.GetSymbol()
        if element == 'Ir':
            target_index = atom.GetIdx()
            break
    assert target_index is not None
    # Then, retrieve the coordinates of the iridium atom in the geometry
    geometry = dative_mol.GetConformer()
    target_point = geometry.GetAtomPosition(target_index)
    # Finally, use these coordinates as the center
    CanonicalizeConformer(geometry, center=target_point)
    return dative_mol

def transfer_conformation(mol, substruct, conformer=0):
    '''Given a molecule, and a second molecule which is a substructure of the
    first, assign coordinates to the substructure based on the matching part of
    the original molecule'''
    match = mol.GetSubstructMatch(substruct)
    substruct_conformation = Conformer(substruct.GetNumAtoms())
    for i, index in enumerate(match):
        point = mol.GetConformer(conformer).GetAtomPosition(index)
        substruct_conformation.SetAtomPosition(i, point)
    substruct.AddConformer(substruct_conformation)

fac = load_template("OHUZEW.mol2")
mer = load_template("OHUZIA.mol2")

template = MolFromSmarts("[Ir]1~n:[*]~[*]:c~1")

carbene_fac = load_template("MAXYIU.mol2")
carbene_mer = load_template("MAXYOA.mol2")

facmer_skeleton_smarts = (
    "[Ir]123"
    "(~[n]~[a]~[a]~[c]~1)"
    "(~[n]~[a]~[a]~[c]~2)"
    "(~[n]~[a]~[a]~[c]~3)"
)

fac_skeleton = MolFromSmarts(facmer_skeleton_smarts)
transfer_conformation(fac, fac_skeleton)

mer_skeleton = MolFromSmarts(facmer_skeleton_smarts)
transfer_conformation(mer, mer_skeleton)

carbene_skeleton_smarts = "[Ir]135(<-[CH0](~N(~*)~*~2)~N(~*~2)~c~c~1)(<-[CH0](~N(~*)~*~4)~N(~*~4)~c~c~3)(<-[CH0](~N(~*)~*~6)~N(~*~6)~c~c~5)"
carbene_fac_skeleton = MolFromSmarts(carbene_skeleton_smarts)
transfer_conformation(carbene_fac, carbene_fac_skeleton)
carbene_mer_skeleton = MolFromSmarts(carbene_skeleton_smarts)
transfer_conformation(carbene_mer, carbene_mer_skeleton)

fac_skeletons = [fac_skeleton, carbene_fac_skeleton]
mer_skeletons = [mer_skeleton, carbene_mer_skeleton]

# Skeletons for tridentate carbenes
# I may have to remake these later if I want to control the isomers
# For now I think it doesn't matter because all the carbene ligands are symmetric?
# For homoleptic:
biplet = load_template("BIPLET.mol2")
biplet_skeleton = MolFromSmarts('[Ir]1234(~[#6](~[#7](~[#6])~[#6])~[#7]~c~c~1~c~[#7]~[#6](~[#7](~[#6])~[#6])~2)~[#6](~[#7](~[#6])~[#6])~[#7]~c~c~3~c~[#7]~[#6](~[#7](~[#6])~[#6])~4')
transfer_conformation(biplet, biplet_skeleton)
# For heteroleptic, three with counterligands of different size
soynom = load_template("SOYNOM.mol2")
soynom_skeleton = MolFromSmarts('[Ir]1234(~*~*~*~*~1~*~*~*~2)~[#6](~[#7](~[#6])~[#6])~[#7]~c~c~3~c~[#7]~[#6](~[#7](~[#6])~[#6])~4')
transfer_conformation(soynom, soynom_skeleton)
uyokur = load_template("UYOKUR.mol2")
uyokur_skeleton = MolFromSmarts('[Ir]1234(~*~*~*~*~*~1~*~*~*~2)~[#6](~[#7](~[#6])~[#6])~[#7]~c~c~3~c~[#7]~[#6](~[#7](~[#6])~[#6])~4')
transfer_conformation(uyokur, uyokur_skeleton)
egufiz = load_template("EGUFIZ.mol2")
egufiz_skeleton = MolFromSmarts('[Ir]1234(~*~*~*~*~*~1~*~*~*~*~2)~[#6](~[#7](~[#6])~[#6])~[#7]~c~c~3~c~[#7]~[#6](~[#7](~[#6])~[#6])~4')
transfer_conformation(egufiz, egufiz_skeleton)
# Homoleptic has to go first, since a later pattern can cover it
tridentate_skeletons = [biplet_skeleton, soynom_skeleton, uyokur_skeleton, egufiz_skeleton]

def octahedral_embed(mol, isomer, clearConfs=True):
    # Make a copy of the molecule to avoid side effects (in particular, removing stereochemistry)
    work = Mol(mol)
    # Needed for some of the mol2 files I got from CSD
    # Will not be able to embed with stereochemistry
    RemoveStereochemistry(work)
    if isomer == "fac":
        skeletons = fac_skeletons
    elif isomer == "mer":
        skeletons = mer_skeletons
    elif isomer == "tridentate":
        skeletons = tridentate_skeletons
    else:
        raise ValueError(f"Isomer should be \"mer\" or \"fac\", given {isomer}")
    finished = False
    for skeleton in skeletons:
        if len(work.GetSubstructMatch(skeleton)) > 0:
            # Carbene embedding with a large template gives output "Could not
            # triangle bounds smooth molecule" and raises a ValueError. But
            # with a small template the imidazole is horribly twisted, probably
            # because it thinks the atoms are aliphatic. Ignoring smoothing
            # failures with the large template, it works
            ConstrainedEmbed(work, skeleton, ignoreSmoothingFailures=True, clearConfs=True)
            finished = True
            break
    if not finished:
        raise ValueError("Doesn't match templates")
    # Copy the conformer back to the original molecule
    new_conf = work.GetConformer()
    if clearConfs:
        mol.RemoveAllConformers()
    return mol.AddConformer(new_conf, assignId=True)
