"""Chemical element data used by the molecule editor.

Colors follow the Jmol/CPK convention, covalent radii are in angstroms
(Cordero et al. 2008). Unknown elements fall back to a gray placeholder.
"""

# symbol: (covalent radius, color)
ELEMENTS = {
    "H":  (0.31, "#f0f0f0"),
    "B":  (0.84, "#ffb5b5"),
    "C":  (0.76, "#a0a4ab"),
    "N":  (0.71, "#4a6cf7"),
    "O":  (0.66, "#ff4d4d"),
    "F":  (0.57, "#90e050"),
    "Na": (1.66, "#ab5cf2"),
    "Mg": (1.41, "#8aff00"),
    "Si": (1.11, "#f0c8a0"),
    "P":  (1.07, "#ff8000"),
    "S":  (1.05, "#ffd633"),
    "Cl": (1.02, "#1ff01f"),
    "K":  (2.03, "#8f40d4"),
    "Ca": (1.76, "#3dff00"),
    "Fe": (1.32, "#e06633"),
    "Cu": (1.32, "#c88033"),
    "Zn": (1.22, "#7d80b0"),
    "Br": (1.20, "#a62929"),
    "Ag": (1.45, "#c0c0c0"),
    "I":  (1.39, "#b040b0"),
    "Au": (1.36, "#ffd123"),
}

FALLBACK = (1.0, "#7a8290")

# elements offered in the editor's quick picker
COMMON = ["C", "H", "O", "N", "S", "F", "Cl", "Br", "P", "B", "Si"]

BOND_TOLERANCE = 1.15  # bond when distance < tolerance * (r1 + r2)


def covalent_radius(symbol: str) -> float:
    return ELEMENTS.get(symbol, FALLBACK)[0]


def element_color(symbol: str) -> str:
    return ELEMENTS.get(symbol, FALLBACK)[1]


def normalize_symbol(text: str) -> str:
    """'cl' -> 'Cl'; returns the input stripped/capitalized, may be unknown."""
    t = text.strip()
    if not t:
        return ""
    return t[0].upper() + t[1:].lower()
