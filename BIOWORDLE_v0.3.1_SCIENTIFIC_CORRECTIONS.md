# BioWordle v0.3.1 — Scientific corrections

Audit date: 2026-11-14. Baseline: `83f05f96102abea5d2fcce0ed75dbf3535127f71`. This evidence records the proposed corrections; it does not certify a release or migration by itself. No additional compounds are proposed.

## Method and limits

The seven `baseline_legacy` input structures and metadata in `data/chemle_compounds.json` were compared against the named PubChem records below; ChEBI is an additional identity/ontology source where cited. PubChem and ChEBI are not necessarily independent primary experimental determinations, and PubChem name/synonym aggregation is not independent proof of structure. RDKit parsed the source SMILES and independently calculated formula, formal charge, canonical isomeric SMILES, and standard InChIKey locally. The target keys below reproduced the cited PubChem keys. Standard InChI can normalize certain resonance and mobile-hydrogen representations; matching it alone cannot prove a chosen tautomer or all stereochemical details. The neutral molecular forms are an explicit game representation, not a claim about dominant protonation in a biological environment. No experimental assay or original isolation structure determination was independently reviewed.

## Findings

| Answer | Baseline → reference standard InChIKey | Reference formula | Finding and action |
|---|---|---|---|
| Artemisinin | `BLUAFEHZUWYNDE-GZEPNZDPSA-N` → `BLUAFEHZUWYNDE-NNWCWBAJSA-N` | C15H22O5 | Same InChI connectivity block but different stereochemical layer. Use naturally occurring (+)-artemisinin with the seven stereocenters of [PubChem CID 68827](https://pubchem.ncbi.nlm.nih.gov/compound/68827), not the legacy stereoisomer CID 138059745. Net neutral lactone/endoperoxide. |
| Cholesterol | `HVYWMOMLDIMFJA-PZIGXULQSA-N` → `HVYWMOMLDIMFJA-DPAQBDIFSA-N` | C27H46O | Same connectivity block, different stereochemistry; use the natural cholest-5-en-3β-ol structure of [PubChem CID 5997](https://pubchem.ncbi.nlm.nih.gov/compound/5997). Net neutral; eight specified carbon stereocenters. |
| Geosmin | `SQFNUAVVZSZPMP-XMACVGRKSA-N` → `JLPUXFOGCDVKGO-TUAOUCFPSA-N` | C12H22O | Different connectivity, despite identical formula. Use natural (−)-(4S,4aS,8aR)-geosmin, a neutral tertiary alcohol, [PubChem CID 29746](https://pubchem.ncbi.nlm.nih.gov/compound/29746); use the specific ChEBI reference CHEBI:46702. Three carbon stereocenters are specified. |
| Prodigiosin | `LCNZNMPZJFZRAP-UHFFFAOYSA-N` → `SZXDNGVQRDTJSD-FOWTUZBSSA-N` | C20H25N3O | Baseline is C19H23N3O, has a short alkyl substituent and differs in ring connectivity; it is not explained by a tautomer. Use the neutral, E-specified [PubChem CID 172897611](https://pubchem.ncbi.nlm.nih.gov/compound/172897611) tripyrrole representation, including the pentyl side chain. |
| Pyocyanin | `HAQCALIJITZJCC-UHFFFAOYSA-N` → `YNCMLFHHXWETLD-UHFFFAOYSA-N` | C13H10N2O | The baseline [CID 14360019](https://pubchem.ncbi.nlm.nih.gov/compound/14360019) is 10-methylphenazin-2(10H)-one, a positional isomer, **not** merely a neutral resonance form of pyocyanin. Replace with [CID 6817](https://pubchem.ncbi.nlm.nih.gov/compound/6817) 5-methylphenazin-1-one. Its net-neutral uncharged keto graph and the net-neutral iminium-olate graph in [ChEBI:8653](https://www.ebi.ac.uk/chebi/CHEBI:8653) share the same standard InChIKey: these are legitimate alternative resonance representations of the correct connectivity. Remove the incorrect 2-one alias and correct the legacy explanation. |
| Quinine | `LOUPRKONTZGTKE-AFHBHXEDSA-N` → `LOUPRKONTZGTKE-WZBLMQSHSA-N` | C20H24N2O2 | Same connectivity, different stereochemical layer. The legacy [CID 94175](https://pubchem.ncbi.nlm.nih.gov/compound/94175) is epiquinidine; use neutral free-base (−)-quinine [CID 3034034](https://pubchem.ncbi.nlm.nih.gov/compound/3034034), with its four specified carbon stereocenters. A neutral free base does not imply unprotonated quinine predominates at every pH. RDKit may also report the bridged tertiary nitrogen as a potential unassigned stereocenter; do not confuse that with missing carbon configuration. |
| Violacein | `PJAOXHSAWDWRMO-UHFFFAOYSA-N` → `SHLJIZCPRXXHHZ-UHFFFAOYSA-N` | C20H13N3O3 | Baseline C28H18N4O3 includes an extra indole and is a different constitution, not a tautomer. [CID 11053](https://pubchem.ncbi.nlm.nih.gov/compound/11053) uses a neutral hydroxy-pyrrole/imino-oxindole graph. Its record's ChEBI description discusses E geometry in a different keto/enol representation; the displayed CID 11053 graph has no specified E/Z bond. Treat its stereochemical specification as a representation limitation; do not infer a resolved E/Z assignment solely from its key. Correct the legacy five-ring structural note. |

## Source disagreements

**Source conflicts requiring resolution before asserting full stereochemical identity.** [ChEBI:82758](https://www.ebi.ac.uk/chebi/CHEBI:82758) identifies prodigiosin with the same formula/connectivity block but the opposite methine stereochemical InChI layer (`SZXDNGVQRDTJSD-VBKFSLOCSA-N`) to PubChem CID 172897611 (`SZXDNGVQRDTJSD-FOWTUZBSSA-N`). The chosen PubChem E representation is an explicitly documented choice, not independent corroboration by ChEBI. An original stereochemical determination and curator review are needed to establish the intended naturally occurring stereoisomer; avoid claiming settled consensus.

[ChEBI:131914](https://www.ebi.ac.uk/chebi/CHEBI:131914) describes an E-configured keto representation and reports `XAPNKXIRQFHCHN-QGOAFFKASA-N`, which differs even in the first block from the PubChem CID 11053 hydroxy-pyrrole representation `SHLJIZCPRXXHHZ-UHFFFAOYSA-N`. The ChEBI page's displayed SMILES includes an indole `n` without an explicitly displayed hydrogen, and the source InChI also has a different hydrogen-layer allocation; the discrepancy is not safely resolved by saying the two are necessarily equivalent tautomers. The PubChem representation is a documented provisional choice for the graph; verify structure, valence and E/Z against authoritative experimental material and curator review before claiming the issue settled.

## Reproducible reference inputs

The following are PubChem isomeric SMILES as retrieved on the audit date. They are reference inputs, not hand-edited canonical strings. RDKit `Chem.MolFromSmiles`, `rdMolDescriptors.CalcMolFormula`, `Chem.GetFormalCharge`, and `Chem.MolToInchiKey` yielded the formula and key in the table and net formal charge zero for each. The project compiler's own canonicalization must determine its exact stored SMILES and tests must assert it.

```text
artemisinin C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C
cholesterol C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C
geosmin C[C@H]1CCC[C@@]2([C@@]1(CCCC2)O)C
prodigiosin CCCCCC1=C/C(=C\C2=C(C=C(N2)C3=CC=CN3)OC)/N=C1C
pyocyanin CN1C2=CC=CC=C2N=C3C1=CC=CC3=O
quinine COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O
violacein C1=CC2=C(C(=O)N=C2C=C1)C3=C(NC(=C3)C4=CNC5=C4C=C(C=C5)O)O
```

The InChI first block is useful for triage, not a complete graph-isomorphism proof. Compare mapped graphs and explicit atom/bond stereochemistry where needed; tautomer or resonance normalization can alter apparent connectivity/bond-order interpretation. These corrections materially change the structures or answer identity, even when the old and new formulas coincide. Historical feedback must remain associated with its original scientific data and rules. A migration must not silently re-score submitted guesses with the new data.
