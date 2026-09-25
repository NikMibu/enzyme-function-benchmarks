"""Sequence-only prompts. The state is the raw amino-acid sequence and its length: no accession,
entry name, protein name, organism, homologs or annotation of any kind reaches a model."""

EC1_CLASSES = {
    "1": ("oxidoreductase", "EC 1: catalyses oxidation-reduction reactions (transfer of electrons, "
                            "hydrogen or oxygen between molecules)"),
    "2": ("transferase", "EC 2: transfers a functional group such as a methyl, acyl, amino, "
                         "phosphate or glycosyl group from one molecule to another"),
    "3": ("hydrolase", "EC 3: cleaves bonds by hydrolysis, e.g. ester, peptide or glycosidic bonds"),
    "4": ("lyase", "EC 4: cleaves C-C, C-O, C-N or other bonds by elimination rather than hydrolysis "
                   "or oxidation, often forming a double bond or ring, or adds groups to double bonds"),
    "5": ("isomerase", "EC 5: rearranges atoms within a single molecule (racemases, epimerases, "
                       "mutases, cis-trans isomerases)"),
    "6": ("ligase", "EC 6: joins two molecules with a new covalent bond, driven by hydrolysis of ATP "
                    "or a similar nucleoside triphosphate"),
}
LABEL_TO_EC1 = {name: cls for cls, (name, _) in EC1_CLASSES.items()}


def state(seq: str) -> str:
    return f"Protein amino-acid sequence ({len(seq)} residues, one-letter code):\n{seq}"


def ec1_question(context: bool = False) -> dict:
    basis = ("Based on this amino-acid sequence and the properties computed from it" if context
             else "Based only on this amino-acid sequence")
    return {"ec1": {
        "type": "choice",
        "instructions": (f"{basis}, which top-level Enzyme Commission (EC) class does this "
                         "enzyme belong to?"),
        "criteria": {name: desc for name, desc in EC1_CLASSES.values()},
    }}


def ec4_question(ecs: list, names: dict) -> dict:
    return {"ec4": {
        "type": "choice",
        "instructions": ("Based only on this amino-acid sequence, which of these enzyme activities "
                         "(Enzyme Commission numbers) does this enzyme have?"),
        "criteria": {f"EC {ec}": names[ec] for ec in sorted(ecs)},
    }}


def label_to_class(stage: str, label: str) -> str:
    return LABEL_TO_EC1[label] if stage == "ec1" else label.removeprefix("EC ")
