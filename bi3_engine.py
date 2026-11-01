"""
bi3_engine.py
=============
Core calculation engine for Bi3+ luminescence-structure relationships in inorganic solids.
Implements the step-by-step nephelauxetic function he(M) / he(Bi3+) and outer-sphere
charge transfer (MMCT & IVCT) models from:
  Philippe Boutinaud, "Luminescence-structure relationships in solids doped with Bi3+",
  Phys. Chem. Chem. Phys., 2023, 25, 11027-11054 (DOI: 10.1039/d3cp00289f).
"""

import json
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd

# ==============================================================================
# 1. EMBEDDED REFERENCE DATABASES (Brese & O'Keeffe BVP, Table 1, Table 2, etc.)
# ==============================================================================

# Matrix III: Valence State (VS) and Valence Number (VN) of anions (p. 11031)
ANION_MATRIX_III: Dict[str, Dict[str, float]] = {
    "F":  {"VS": 1.0, "VN": 7.0, "group": "VII"},
    "Cl": {"VS": 1.0, "VN": 7.0, "group": "VII"},
    "Br": {"VS": 1.0, "VN": 7.0, "group": "VII"},
    "I":  {"VS": 1.0, "VN": 7.0, "group": "VII"},
    "O":  {"VS": 2.0, "VN": 6.0, "group": "VI"},
    "S":  {"VS": 2.0, "VN": 6.0, "group": "VI"},
    "Se": {"VS": 2.0, "VN": 6.0, "group": "VI"},
    "Te": {"VS": 2.0, "VN": 6.0, "group": "VI"},
    "N":  {"VS": 3.0, "VN": 5.0, "group": "V"},
}

# Brese & O'Keeffe (1991) Bond Valence Parameters R_ij (in Angstroms)
# Keyed by (Cation_Element, Anion_Element)
BVP_DATABASE: Dict[Tuple[str, str], float] = {
    # Oxides
    ("Li", "O"): 1.466, ("Na", "O"): 1.803, ("K", "O"):  2.132, ("Rb", "O"): 2.263, ("Cs", "O"): 2.417,
    ("Mg", "O"): 1.693, ("Ca", "O"): 1.967, ("Sr", "O"): 2.118, ("Ba", "O"): 2.285, ("Zn", "O"): 1.704,
    ("Cd", "O"): 1.904, ("Pb", "O"): 2.112, ("B", "O"):  1.371, ("Al", "O"): 1.651, ("Ga", "O"): 1.730,
    ("In", "O"): 1.902, ("Sc", "O"): 1.849, ("Y", "O"):  2.014, ("La", "O"): 2.172, ("Gd", "O"): 2.065,
    ("Lu", "O"): 1.971, ("Bi", "O"): 2.094, ("Si", "O"): 1.624, ("Ge", "O"): 1.748, ("Sn", "O"): 1.905,
    ("Ti", "O"): 1.815, ("Zr", "O"): 1.928, ("Hf", "O"): 1.923, ("P", "O"):  1.617, ("V", "O"):  1.803,
    ("As", "O"): 1.767, ("Nb", "O"): 1.911, ("Sb", "O"): 1.973, ("Ta", "O"): 1.920, ("S", "O"):  1.624,
    ("Mo", "O"): 1.907, ("W", "O"):  1.917,
    # Fluorides
    ("Li", "F"): 1.360, ("Na", "F"): 1.677, ("K", "F"):  1.992, ("Rb", "F"): 2.160, ("Cs", "F"): 2.330,
    ("Mg", "F"): 1.581, ("Ca", "F"): 1.842, ("Sr", "F"): 2.019, ("Ba", "F"): 2.188, ("Y", "F"):  1.904,
    ("La", "F"): 2.057, ("Gd", "F"): 1.950, ("Lu", "F"): 1.860, ("Sc", "F"): 1.731, ("Al", "F"): 1.545,
    # Chlorides
    ("Na", "Cl"): 2.150, ("K", "Cl"): 2.520, ("Rb", "Cl"): 2.650, ("Cs", "Cl"): 2.790,
    ("Mg", "Cl"): 2.080, ("Ca", "Cl"): 2.370, ("Sr", "Cl"): 2.510, ("Ba", "Cl"): 2.690,
    ("Y", "Cl"):  2.400, ("La", "Cl"): 2.545, ("Gd", "Cl"): 2.445, ("Zr", "Cl"): 2.330,
    ("Hf", "Cl"): 2.320, ("Sn", "Cl"): 2.360, ("Bi", "Cl"): 2.480,
    # Bromides & Iodides
    ("K", "Br"):  2.660, ("Cs", "Br"): 2.950, ("Ca", "Br"): 2.490, ("Cd", "Br"): 2.450,
    ("Y", "Br"):  2.550, ("La", "Br"): 2.720, ("Ca", "I"):  2.720, ("K", "I"):   2.880,
    # Chalcogenides
    ("Mg", "S"):  2.130, ("Ca", "S"):  2.450, ("Sr", "S"):  2.590, ("Ba", "S"):  2.770,
    ("Zn", "S"):  2.170, ("La", "S"):  2.640, ("Ca", "Se"): 2.560, ("Ca", "Te"): 2.760,
}

# Table 1 of Boutinaud (2023): Local distortion correction Delta_d (in Angstroms) after Bi3+ substitution
TABLE1_DELTA_D: Dict[str, float] = {
    "Cs": -0.285, "Rb": -0.245, "K":  -0.175, "Ba": -0.125, "Pb": -0.070,
    "Sr": -0.050, "Na":  0.000, "Ca":  0.000, "La":  0.000, "Cd": +0.035,
    "Gd": +0.060, "Eu": +0.055, "Sm": +0.045, "Nd": +0.025, "Y":  +0.075,
    "Lu": +0.090, "In": +0.120, "Li": +0.135, "Zn": +0.140, "Mg": +0.150,
    "Sc": +0.150, "Sn": +0.175, "Ga": +0.200, "Zr": +0.155, "Hf": +0.160,
    "Ti": +0.210, "Al": +0.245, "Bi":  0.000,
}

# Table 2 of Boutinaud (2023, PDF p. 4): "Value of parameters k_CN^CN' and a_CN^CN' in eqn (4) for MMCT and IVCT
# transitions as a function of CN and CN'" (columns: chi_CN(Bi3+); MMCT with CN' = 4: k, a; MMCT with CN' >= 6: k, a;
# IVCT: k, a). The paper lists CN = 5, 6, 7, 8, 9, 10 and 12 only.
# The row for CN = 11 is NOT in the paper: it is the linear interpolation (midpoints) of the CN = 10 and CN = 12 rows.
# Sites with CN > 12 use the CN = 12 row (and CN < 5 the CN = 5 row); an acceptor with CN' other than 4 uses the
# CN' >= 6 columns. The application flags such sites (bi3_export.ct_parameter_note).
TABLE2_CT_PARAMS: Dict[int, Dict[str, float]] = {
    5:  {"chi_Bi": 1.430, "k_CN4": 6.07, "a_CN4": 1.06, "k_CN6": 4.77, "a_CN6": 1.18, "k_IVCT": 8.00, "a_IVCT": 2.58},
    6:  {"chi_Bi": 1.400, "k_CN4": 6.20, "a_CN4": 1.04, "k_CN6": 4.87, "a_CN6": 1.16, "k_IVCT": 8.16, "a_IVCT": 2.52},
    7:  {"chi_Bi": 1.370, "k_CN4": 6.33, "a_CN4": 1.02, "k_CN6": 4.98, "a_CN6": 1.13, "k_IVCT": 8.34, "a_IVCT": 2.46},
    8:  {"chi_Bi": 1.340, "k_CN4": 6.48, "a_CN4": 0.99, "k_CN6": 5.09, "a_CN6": 1.11, "k_IVCT": 8.53, "a_IVCT": 2.41},
    9:  {"chi_Bi": 1.330, "k_CN4": 6.52, "a_CN4": 0.99, "k_CN6": 5.13, "a_CN6": 1.10, "k_IVCT": 8.59, "a_IVCT": 2.39},
    10: {"chi_Bi": 1.305, "k_CN4": 6.65, "a_CN4": 0.97, "k_CN6": 5.22, "a_CN6": 1.08, "k_IVCT": 8.76, "a_IVCT": 2.35},
    11: {"chi_Bi": 1.292, "k_CN4": 6.71, "a_CN4": 0.96, "k_CN6": 5.27, "a_CN6": 1.07, "k_IVCT": 8.84, "a_IVCT": 2.32},
    12: {"chi_Bi": 1.280, "k_CN4": 6.78, "a_CN4": 0.95, "k_CN6": 5.33, "a_CN6": 1.06, "k_IVCT": 8.93, "a_IVCT": 2.30},
}

# Figure 2 of Boutinaud (2023) & Li & Xue (2006, J. Phys. Chem. A 110, 11332): electronegativities chi_CN'(acceptor)
# (the 2023 model uses plain electronegativities; optical ones belong to the 2011/2013 models)
ACCEPTOR_CHI_DATABASE: Dict[Tuple[str, int], float] = {
    # CN' >= 6 (Octahedral / High CN)
    ("Sc", 6): 1.41, ("Y", 6):  1.34, ("La", 6): 1.33, ("Lu", 6): 1.37, ("Gd", 6): 1.35,
    ("Y", 7):  1.31, ("Y", 8):  1.28, ("Y", 9):  1.34, ("La", 8): 1.25, ("La", 9): 1.24, ("La", 12): 1.20,
    ("Lu", 7): 1.34, ("Lu", 8): 1.31, ("Sc", 8): 1.35,
    ("Ti", 6): 1.73, ("Zr", 6): 1.71, ("Hf", 6): 1.61,
    ("V", 6):  2.03, ("Nb", 6): 1.92, ("Ta", 6): 1.86,
    ("Mo", 6): 2.17, ("W", 6):  2.10,
    ("Cu", 6): 1.16, ("Ag", 6): 1.09, ("Zn", 6): 1.33, ("Cd", 6): 1.28,
    ("Ga", 6): 1.58, ("In", 6): 1.48, ("Ge", 6): 1.85, ("Sn", 6): 1.71,
    ("As", 6): 2.16, ("Sb", 6): 1.97,
    # CN' = 4 (Tetrahedral)
    ("V", 4):  2.46, ("Nb", 4): 2.12, ("Mo", 4): 2.59, ("W", 4):  2.50,
    ("Cu", 4): 1.23, ("Ag", 4): 1.12, ("Zn", 4): 1.35, ("Cd", 4): 1.35,
    ("Ga", 4): 1.75, ("In", 4): 1.63, ("Ge", 4): 2.12, ("Sn", 4): 1.87,
    ("As", 4): 2.49,
}

ATOMIC_MASSES: Dict[str, float] = {
    "H": 1.008, "Li": 6.94, "B": 10.81, "N": 14.007, "O": 15.999, "F": 18.998,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Rb": 85.468,
    "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Ag": 107.87,
    "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Eu": 151.96, "Gd": 157.25, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Pb": 207.2, "Bi": 208.98, "Zn": 65.38,
}


def strip_element_symbol(label: str) -> str:
    """Extracts the clean chemical element symbol from a site label (e.g., 'Ca1' -> 'Ca', 'Y2+' -> 'Y')."""
    m = re.match(r"^\s*([A-Z][a-z]?)", str(label))
    return m.group(1) if m else str(label).strip()


def get_bvp(cation_label: str, anion_label: str) -> float:
    cat_el = strip_element_symbol(cation_label)
    an_el = strip_element_symbol(anion_label)
    return BVP_DATABASE.get((cat_el, an_el), 1.950)


def get_delta_d(cation_label: str) -> float:
    cat_el = strip_element_symbol(cation_label)
    return TABLE1_DELTA_D.get(cat_el, 0.0)


def get_acceptor_chi(acceptor_label: str, cn_prime: int) -> Optional[float]:
    el = strip_element_symbol(acceptor_label)
    if (el, cn_prime) in ACCEPTOR_CHI_DATABASE:
        return ACCEPTOR_CHI_DATABASE[(el, cn_prime)]
    if cn_prime >= 6 and (el, 6) in ACCEPTOR_CHI_DATABASE:
        return ACCEPTOR_CHI_DATABASE[(el, 6)]
    if cn_prime == 4 and (el, 4) in ACCEPTOR_CHI_DATABASE:
        return ACCEPTOR_CHI_DATABASE[(el, 4)]
    return None


# ==============================================================================
# 2. CORE BOUTINAUD (2023) MATHEMATICAL ENGINE (Equations 1 - 19)
# ==============================================================================

@dataclass
class CrystalStructureInput:
    name: str
    formula: str
    space_group: str
    Z: int
    V: float  # Unit cell volume in A^3
    # sites_cations: site_label -> {'element': str, 'w': float, 'N': int, 'BVP': {an_label: float}, 'Q_override': Optional[float], 'bonds': {an_label: {'N_MX': int, 'd_list': List[float]}}}
    sites_cations: Dict[str, Dict[str, Any]]
    # sites_anions: an_label -> {'element': str, 'w': float, 'N': int, 'bonds': {cat_label: int}}
    sites_anions: Dict[str, Dict[str, Any]]
    # acceptors: list of dicts {'doping_site': str, 'acceptor_site': str, 'acceptor_element': str, 'CN_prime': int, 'd_host': float}
    acceptors: List[Dict[str, Any]] = field(default_factory=list)
    # optional experimental data from archival literature or user input
    exp_EA_eV: Optional[float] = None
    exp_S_eV: Optional[float] = None
    pub_he_M: Optional[Dict[str, float]] = None
    pub_he_Bi: Optional[Dict[str, float]] = None
    reference_note: str = ""


def run_boutinaud_2023_calculation(
    crystal: CrystalStructureInput,
    custom_exp_ex_nm: Optional[float] = None,
    custom_exp_em_nm: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Executes the complete step-by-step calculation of he(M), he(Bi3+), EA, EC, EMMCT, EIVCT,
    and Delta_E(Bi3+) for all candidate doping sites in `crystal`.
    """
    Z = float(crystal.Z)
    V = float(crystal.V)

    # Step 1: Compute Bond Valence Sum Q(M) for every cation site (Eq. 10)
    cation_Q: Dict[str, float] = {}
    for c_label, c_data in crystal.sites_cations.items():
        if c_data.get("Q_override") is not None:
            cation_Q[c_label] = float(c_data["Q_override"])
        else:
            q_val = 0.0
            for x_label, b_info in c_data["bonds"].items():
                bvp = c_data.get("BVP", {}).get(x_label, get_bvp(c_data["element"], crystal.sites_anions[x_label]["element"]))
                d_arr = np.array(b_info["d_list"], dtype=float)
                q_val += float(np.sum(np.exp(-(d_arr - bvp) / 0.37)))
            cation_Q[c_label] = q_val

    # Build Matrix I & Matrix II summary rows
    matrix_I_rows = []
    for c_label, c_data in crystal.sites_cations.items():
        row = {
            "Site": c_label,
            "Role": "Cation",
            "Element": c_data["element"],
            "Multiplicity w": c_data["w"],
            "Stoich. Coeff (w/Z)": round(c_data["w"] / Z, 4),
            "Total CN (N)": c_data["N"],
            "BVS Q(M) [Eq. 10]": round(cation_Q[c_label], 4),
            "Delta_d (A) [Table 1]": get_delta_d(c_data["element"]),
        }
        for x_label in crystal.sites_anions:
            n_mx = c_data["bonds"].get(x_label, {}).get("N_MX", 0)
            row[f"N(Site-{x_label})"] = n_mx
        matrix_I_rows.append(row)

    for x_label, x_data in crystal.sites_anions.items():
        row = {
            "Site": x_label,
            "Role": "Anion",
            "Element": x_data["element"],
            "Multiplicity w": x_data["w"],
            "Stoich. Coeff (w/Z)": round(x_data["w"] / Z, 4),
            "Total CN (N)": x_data["N"],
            "BVS Q(M) [Eq. 10]": np.nan,
            "Delta_d (A) [Table 1]": np.nan,
        }
        for c_label in crystal.sites_cations:
            n_xm = x_data["bonds"].get(c_label, 0)
            row[f"N(Site-{c_label})"] = n_xm
        matrix_I_rows.append(row)

    # Helper to compute denominator of Eq. 16 for either unperturbed host (delta_d = 0)
    # or with local distortion delta_d at active_site
    def compute_eq16_denom(active_site: Optional[str], delta_d_val: float) -> float:
        denom = 0.0
        for c_label, c_data in crystal.sites_cations.items():
            w_C = float(c_data["w"])
            dd = delta_d_val if c_label == active_site else 0.0
            for x_label, b_info in c_data["bonds"].items():
                N_CX = float(b_info["N_MX"])
                if N_CX > 0:
                    d_mean = float(np.mean(b_info["d_list"])) + dd
                    denom += (w_C * N_CX / V) * (d_mean ** 3)
        return denom

    denom_host = compute_eq16_denom(None, 0.0)

    # Step 2: Evaluate sub-shell dielectric parameters (Eqs. 7 - 17) for each cation site
    subshell_rows = []
    site_summary_rows = []

    exp_EA_used = (1239.84193 / custom_exp_ex_nm) if (custom_exp_ex_nm and custom_exp_ex_nm > 0) else crystal.exp_EA_eV
    exp_Em_used = (1239.84193 / custom_exp_em_nm) if (custom_exp_em_nm and custom_exp_em_nm > 0) else (
        (crystal.exp_EA_eV - crystal.exp_S_eV) if (crystal.exp_EA_eV and crystal.exp_S_eV) else None
    )
    exp_S_used = (exp_EA_used - exp_Em_used) if (exp_EA_used is not None and exp_Em_used is not None) else crystal.exp_S_eV

    for c_label, c_data in crystal.sites_cations.items():
        N_M = int(c_data["N"])
        w_M = float(c_data["w"])
        Q_M = cation_Q[c_label]
        delta_d = float(c_data.get("delta_d_override", get_delta_d(c_data["element"])))
        denom_doped = compute_eq16_denom(c_label, delta_d)

        he_M_sq = 0.0
        he_Bi_sq = 0.0
        all_host_d = []

        for x_label, b_info in c_data["bonds"].items():
            N_MX = float(b_info["N_MX"])
            if N_MX <= 0:
                continue
            d_host = float(np.mean(b_info["d_list"]))
            all_host_d.extend(b_info["d_list"])
            d_Bi = d_host + delta_d

            x_data = crystal.sites_anions[x_label]
            w_X = float(x_data["w"])
            N_X = float(x_data["N"])
            N_XM = float(x_data["bonds"][c_label])
            el_X = strip_element_symbol(x_data["element"])
            VS_X = ANION_MATRIX_III.get(el_X, {"VS": 2.0})["VS"]
            VN_X = ANION_MATRIX_III.get(el_X, {"VN": 6.0})["VN"]

            # Kappa functions (Eq. 7)
            K_MX = (w_M / Z) * (N_MX / N_M)
            K_XM = (w_X / Z) * (N_XM / N_X)

            # Effective anion charge Q(Xi) (Eq. 9)
            Q_Xi = (K_MX / K_XM) * Q_M

            # Gamma_2 (Eq. 13) & N* (Eq. 15)
            Gamma_2 = Q_M * (1.0 - VN_X / VS_X)
            N_star = (Q_M / N_M) + (VN_X * K_MX * Q_M) / (VS_X * K_XM * N_X)

            # Evaluate both Unrelaxed Host (d_host, denom_host) and Bi3+-Doped (d_Bi, denom_doped)
            num_G1_base = (K_MX * N_M + K_XM * N_X) / (K_MX + K_XM)

            # Host state (M)
            G1_M = 2.56 * (num_G1_base ** 2) / d_host
            v_M = (d_host ** 3) / denom_host
            G3_M = math.exp(-1.364 * d_host * ((N_star / v_M) ** (1.0 / 6.0)))
            fc_M = 1.0 / (1.0 + ((G1_M * Gamma_2 * G3_M) ** 2) / (1579.0 * (d_host ** -4.96)))
            al_M = (0.208 * fc_M * (d_host ** 4.96) * N_star) / (
                3.0 + 0.874 * fc_M * (d_host ** 4.96) * (N_star / v_M)
            )
            he_M_sq += N_MX * fc_M * al_M * (Q_Xi ** 2)

            # Bi3+-doped state (Bi3+)
            G1_Bi = 2.56 * (num_G1_base ** 2) / d_Bi
            v_Bi = (d_Bi ** 3) / denom_doped
            G3_Bi = math.exp(-1.364 * d_Bi * ((N_star / v_Bi) ** (1.0 / 6.0)))
            fc_Bi = 1.0 / (1.0 + ((G1_Bi * Gamma_2 * G3_Bi) ** 2) / (1579.0 * (d_Bi ** -4.96)))
            al_Bi = (0.208 * fc_Bi * (d_Bi ** 4.96) * N_star) / (
                3.0 + 0.874 * fc_Bi * (d_Bi ** 4.96) * (N_star / v_Bi)
            )
            he_Bi_sq += N_MX * fc_Bi * al_Bi * (Q_Xi ** 2)

            subshell_rows.append({
                "Cation Site": c_label,
                "Anion Site": x_label,
                "N(M-Xi)": int(N_MX),
                "N(Xi-M)": int(N_XM),
                "d_host (A)": round(d_host, 4),
                "d_Bi (A)": round(d_Bi, 4),
                "K(M-Xi) [Eq. 7]": round(K_MX, 4),
                "K(Xi-M) [Eq. 7]": round(K_XM, 4),
                "Q(Xi) [Eq. 9]": round(Q_Xi, 4),
                "N*(M-Xi) [Eq. 15]": round(N_star, 4),
                "v_host (A^3) [Eq. 16]": round(v_M, 4),
                "Gamma_1 (host) [Eq. 12]": round(G1_M, 3),
                "Gamma_2 [Eq. 13]": round(Gamma_2, 3),
                "Gamma_3 (host) [Eq. 14]": round(G3_M, 5),
                "f_c (host) [Eq. 11]": round(fc_M, 5),
                "alpha (host, A^3) [Eq. 17]": round(al_M, 4),
                "f_c (Bi3+)": round(fc_Bi, 5),
                "alpha (Bi3+, A^3)": round(al_Bi, 4),
            })

        he_M = math.sqrt(he_M_sq)
        he_Bi = math.sqrt(he_Bi_sq)

        # Only sites with CN >= 5 can accommodate Bi3+ in Boutinaud's model (Table 2)
        is_viable_bi_site = (N_M >= 5) and (strip_element_symbol(c_data["element"]) not in ("P", "B", "Si", "S", "V", "Mo", "W", "As"))

        # Step 3: Intra-ionic A-band (Eq. 18) and C-band (Eq. 19)
        EA_eV = 9.41 * (1.0 - 0.555 * (he_Bi ** 0.240)) if he_Bi > 0 else 9.41
        EA_nm = 1239.84193 / EA_eV
        EA_cm = EA_eV * 8065.54

        EC_eV = 3.236 + 2.29 * ((max(EA_eV - 2.972, 0.0)) ** 0.856)
        EC_nm = 1239.84193 / EC_eV

        # Step 4: Extra-ionic MMCT and Bi3+-Bi3+ Pair IVCT (Eq. 4 & Table 2)
        cn_key = min(max(N_M, 5), 12)
        if cn_key not in TABLE2_CT_PARAMS:
            cn_key = 9
        t2 = TABLE2_CT_PARAMS[cn_key]
        chi_Bi = t2["chi_Bi"]

        # Find lowest MMCT energy among all acceptors linked to this doping site
        site_acceptors = [acc for acc in crystal.acceptors if acc["doping_site"] == c_label]
        best_mmct_eV = None
        best_mmct_label = "None"
        best_ratio_chi_d = None
        best_ivct_eV = None
        best_ivct_d = None

        for acc in site_acceptors:
            cn_p = int(acc["CN_prime"])
            d_host_acc = float(acc["d_host"])
            d_bi_acc = d_host_acc + delta_d
            acc_el = strip_element_symbol(acc["acceptor_element"])

            # Check if acceptor is another Bi3+ ion (IVCT pair) or a host d0/d10 cation (MMCT)
            if acc_el == "Bi":
                e_ivct = t2["k_IVCT"] * (chi_Bi - t2["a_IVCT"] * (chi_Bi / d_bi_acc))
                if best_ivct_eV is None or e_ivct < best_ivct_eV:
                    best_ivct_eV = e_ivct
                    best_ivct_d = d_bi_acc
            else:
                chi_acc = acc.get("chi_override") or get_acceptor_chi(acc_el, cn_p)
                if chi_acc is not None:
                    k_val = t2["k_CN4"] if cn_p == 4 else t2["k_CN6"]
                    a_val = t2["a_CN4"] if cn_p == 4 else t2["a_CN6"]
                    e_mmct = k_val * (chi_Bi - a_val * (chi_acc / d_bi_acc))
                    if best_mmct_eV is None or e_mmct < best_mmct_eV:
                        best_mmct_eV = e_mmct
                        best_mmct_label = f"{acc['acceptor_site']}[{cn_p}] (d_host={d_host_acc:.3f} A)"
                        best_ratio_chi_d = chi_acc / d_bi_acc

                    # Also compute Bi-Bi IVCT pair at the same self-cation distance if same sublattice
                    if strip_element_symbol(acc_el) == strip_element_symbol(c_data["element"]):
                        e_ivct_self = t2["k_IVCT"] * (chi_Bi - t2["a_IVCT"] * (chi_Bi / d_bi_acc))
                        if best_ivct_eV is None or e_ivct_self < best_ivct_eV:
                            best_ivct_eV = e_ivct_self
                            best_ivct_d = d_bi_acc

        # Classify Situation (a), (b), or (c) from Fig. 3 of Boutinaud (2023)
        if best_mmct_eV is None:
            situation_str = "Situation (a) [↓]: Pure Intra-Ionic (No d0/d10 MMCT)"
        elif EA_eV < best_mmct_eV - 0.35:
            situation_str = "Situation (a) [↓]: E_A < E_CT - 0.35 eV (A-State Below CT)"
        elif EA_eV > best_mmct_eV + 0.35:
            situation_str = "Situation (c) [↑]: E_A > E_CT + 0.35 eV (CT-State Below A)"
        else:
            situation_str = "Situation (b) [≈]: |E_A - E_CT| <= 0.35 eV (Resonant A & CT)"

        # Step 5: Emission Origin Criterion Delta_E(Bi3+) = |E_CT - E_A - S| (Eq. 5)
        delta_E_crit = None
        emission_origin = "A-State Emission (Intra-Ionic 3P1 -> 1S0)"
        if best_mmct_eV is not None and exp_S_used is not None:
            delta_E_crit = abs(best_mmct_eV - EA_eV - exp_S_used)
            if delta_E_crit < 1.39:
                emission_origin = f"A-State Emission (Delta_E = {delta_E_crit:.2f} eV < 1.39 eV)"
            elif delta_E_crit > 1.41:
                emission_origin = f"CT / Trapped Exciton Emission (Delta_E = {delta_E_crit:.2f} eV > 1.41 eV)"
            else:
                emission_origin = f"Borderline A / CT Regime (Delta_E = {delta_E_crit:.2f} eV in 1.39-1.41 eV)"

        site_summary_rows.append({
            "Doping Site": f"{c_label}[{N_M}]",
            "Raw Site": c_label,
            "Element": c_data["element"],
            "CN": N_M,
            "Viable Bi3+ Site": is_viable_bi_site,
            "Mean d_host (A)": round(float(np.mean(all_host_d)), 4) if all_host_d else np.nan,
            "Delta_d (A)": round(delta_d, 4),
            "Q(M) [BVS]": round(Q_M, 3),
            "he(M) [Calc]": round(he_M, 3),
            "he(M) [Pub]": crystal.pub_he_M.get(c_label) if crystal.pub_he_M else None,
            "he(Bi3+) [Calc]": round(he_Bi, 3),
            "he(Bi3+) [Pub]": crystal.pub_he_Bi.get(c_label) if crystal.pub_he_Bi else None,
            "E_A (eV) [Eq. 18]": round(EA_eV, 2),
            "lambda_A (nm)": round(EA_nm, 1),
            "E_A (cm^-1)": int(round(EA_cm, 0)),
            "E_C (eV) [Eq. 19]": round(EC_eV, 2),
            "lambda_C (nm)": round(EC_nm, 1),
            "Acceptor [CN']": best_mmct_label,
            "chi_CN'/d (A^-1)": round(best_ratio_chi_d, 3) if best_ratio_chi_d else None,
            "E_MMCT (eV) [Eq. 4]": round(best_mmct_eV, 2) if best_mmct_eV else None,
            "lambda_MMCT (nm)": round(1239.84193 / best_mmct_eV, 1) if best_mmct_eV else None,
            "E_IVCT Pair (eV) [Eq. 4]": round(best_ivct_eV, 2) if best_ivct_eV else None,
            "lambda_IVCT Pair (nm)": round(1239.84193 / best_ivct_eV, 1) if best_ivct_eV else None,
            "Exp. E_exc (eV)": round(exp_EA_used, 2) if exp_EA_used else None,
            "Exp. Stokes Shift S (eV)": round(exp_S_used, 2) if exp_S_used else None,
            "Delta_E(Bi3+) (eV) [Eq. 5]": round(delta_E_crit, 2) if delta_E_crit is not None else None,
            "Excited-State Regime": situation_str,
            "Predicted Emission Origin": emission_origin,
        })

    return {
        "crystal": crystal,
        "matrix_I_df": pd.DataFrame(matrix_I_rows),
        "subshells_df": pd.DataFrame(subshell_rows),
        "summary_df": pd.DataFrame(site_summary_rows),
        "denom_host": denom_host,
    }


# ==============================================================================
# 3. CIF AND JSON PARSERS (WITH MISSING-PARAMETER DETECTION FOR JSON)
# ==============================================================================

def _parse_cif_float(val_str: str) -> float:
    """Strips standard uncertainty parentheses e.g. '3.6451(1)' -> 3.6451 or '1/2' -> 0.5."""
    s = str(val_str).strip().strip("'").strip('"')
    s = re.sub(r"\(\d+\)", "", s)
    if "/" in s:
        parts = s.split("/")
        return float(parts[0]) / float(parts[1])
    return float(s)


def _eval_symop_coord(expr: str, x: float, y: float, z: float) -> float:
    """Evaluates a single symmetry coordinate expression like '-x+1/2' or 'y+0.5'."""
    e = expr.strip().lower().replace(" ", "")
    # Replace fractions p/q with float
    e = re.sub(r"(\d+)/(\d+)", lambda m: str(float(m.group(1)) / float(m.group(2))), e)
    val = 0.0
    # Match terms like +x, -x, +0.5, etc.
    for token in re.findall(r"[+-]?[^+-]+", e):
        if "x" in token:
            coeff = token.replace("x", "")
            c = -1.0 if coeff == "-" else (1.0 if coeff in ("", "+") else float(coeff))
            val += c * x
        elif "y" in token:
            coeff = token.replace("y", "")
            c = -1.0 if coeff == "-" else (1.0 if coeff in ("", "+") else float(coeff))
            val += c * y
        elif "z" in token:
            coeff = token.replace("z", "")
            c = -1.0 if coeff == "-" else (1.0 if coeff in ("", "+") else float(coeff))
            val += c * z
        else:
            val += float(token)
    return val % 1.0


def build_crystal_from_unit_cell_sites(
    name: str,
    formula: str,
    space_group: str,
    a: float,
    b: float,
    c: float,
    alpha_deg: float,
    beta_deg: float,
    gamma_deg: float,
    Z: int,
    asymmetric_sites: List[Dict[str, Any]],
    symops: Optional[List[str]] = None,
    max_cutoff_factor: float = 1.32,
) -> CrystalStructureInput:
    """
    Generates all symmetry-equivalent positions in the unit cell, computes the metric tensor,
    builds coordination shells (Matrix I & Matrix II), and identifies shortest cation-acceptor distances.
    """
    # Compute unit cell metric tensor and Cartesian transformation matrix
    ar, br, gr = math.radians(alpha_deg), math.radians(beta_deg), math.radians(gamma_deg)
    cos_a, cos_b, cos_g = math.cos(ar), math.cos(br), math.cos(gr)
    sin_g = math.sin(gr)
    vol_factor = math.sqrt(max(1.0 - cos_a**2 - cos_b**2 - cos_g**2 + 2.0 * cos_a * cos_b * cos_g, 1e-12))
    V = a * b * c * vol_factor

    # Fractional to Cartesian matrix M_cart (columns are lattice vectors a, b, c)
    M_cart = np.array([
        [a, b * cos_g, c * cos_b],
        [0.0, b * sin_g, c * (cos_a - cos_b * cos_g) / sin_g],
        [0.0, 0.0, c * vol_factor / sin_g],
    ])

    # Default symmetry operations if none supplied (handles I4mm / I4/mmm body-centered tetragonal or P1)
    if not symops:
        sg_clean = space_group.replace(" ", "").strip()
        if "I4/mmm" in sg_clean or "139" in sg_clean:
            base_ops = [
                "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z",
                "-x,y,-z", "x,-y,-z", "y,x,-z", "-y,-x,-z",
                "-x,-y,-z", "x,y,-z", "y,-x,-z", "-y,x,-z",
                "x,-y,z", "-x,y,z", "-y,-x,z", "y,x,z"
            ]
            symops = base_ops + [
                f"{op.split(',')[0]}+1/2,{op.split(',')[1]}+1/2,{op.split(',')[2]}+1/2" for op in base_ops
            ]
        elif "I4mm" in sg_clean or "107" in sg_clean:
            base_ops = [
                "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z",
                "x,-y,z", "-x,y,z", "-y,-x,z", "y,x,z"
            ]
            symops = base_ops + [
                f"{op.split(',')[0]}+1/2,{op.split(',')[1]}+1/2,{op.split(',')[2]}+1/2" for op in base_ops
            ]
        else:
            symops = ["x,y,z"]

    # Expand asymmetric sites into full unit cell orbit positions
    expanded_orbits: Dict[str, Dict[str, Any]] = {}
    for idx, site in enumerate(asymmetric_sites):
        el = strip_element_symbol(site["element"])
        label = site.get("label") or f"{el}{idx+1}"
        fx, fy, fz = _parse_cif_float(site["x"]), _parse_cif_float(site["y"]), _parse_cif_float(site["z"])
        orbit_coords = []
        for op in symops:
            parts = op.split(",")
            if len(parts) == 3:
                nx = _eval_symop_coord(parts[0], fx, fy, fz)
                ny = _eval_symop_coord(parts[1], fx, fy, fz)
                nz = _eval_symop_coord(parts[2], fx, fy, fz)
                cand = np.array([nx, ny, nz])
                # Check duplicate modulo 1.0
                is_dup = False
                for existing in orbit_coords:
                    diff = (cand - existing + 0.5) % 1.0 - 0.5
                    if np.linalg.norm(M_cart @ diff) < 0.05:
                        is_dup = True
                        break
                if not is_dup:
                    orbit_coords.append(cand)
        expanded_orbits[label] = {
            "element": el,
            "is_anion": el in ANION_MATRIX_III,
            "w": float(len(orbit_coords)),
            "coords": orbit_coords,
        }

    # Separate cation orbits and anion orbits
    cat_labels = [k for k, v in expanded_orbits.items() if not v["is_anion"]]
    an_labels = [k for k, v in expanded_orbits.items() if v["is_anion"]]

    # Build 3x3x3 supercell translations for periodic minimum-image distance search
    shifts = np.array([[ix, iy, iz] for ix in (-1, 0, 1) for iy in (-1, 0, 1) for iz in (-1, 0, 1)], dtype=float)

    sites_cations: Dict[str, Dict[str, Any]] = {}
    sites_anions: Dict[str, Dict[str, Any]] = {
        x_lab: {"element": expanded_orbits[x_lab]["element"], "w": expanded_orbits[x_lab]["w"], "N": 0, "bonds": {}}
        for x_lab in an_labels
    }

    for c_lab in cat_labels:
        c_info = expanded_orbits[c_lab]
        c_el = c_info["element"]
        c_pos0 = c_info["coords"][0]  # Representative position of Wyckoff orbit

        # Collect all distances to all anion orbits within 3.6 A
        all_an_neighbors = []
        for x_lab in an_labels:
            x_info = expanded_orbits[x_lab]
            for x_pos in x_info["coords"]:
                for sh in shifts:
                    d_vec = (x_pos + sh) - c_pos0
                    dist = float(np.linalg.norm(M_cart @ d_vec))
                    if 0.8 < dist < 3.6:
                        all_an_neighbors.append((dist, x_lab))

        all_an_neighbors.sort(key=lambda t: t[0])
        if not all_an_neighbors:
            continue
        d_min = all_an_neighbors[0][0]
        # Coordination shell cutoff: max_cutoff_factor * d_min (typically 1.32 * d_min captures first polyhedron)
        d_cutoff = max(d_min * max_cutoff_factor, d_min + 0.45)
        coord_neighbors = [(d, xl) for (d, xl) in all_an_neighbors if d <= d_cutoff]

        bonds_dict: Dict[str, Dict[str, Any]] = {}
        bvp_dict: Dict[str, float] = {}
        for dist, x_lab in coord_neighbors:
            if x_lab not in bonds_dict:
                bonds_dict[x_lab] = {"N_MX": 0, "d_list": []}
                bvp_dict[x_lab] = get_bvp(c_el, expanded_orbits[x_lab]["element"])
            bonds_dict[x_lab]["N_MX"] += 1
            bonds_dict[x_lab]["d_list"].append(dist)

        total_cn = sum(b["N_MX"] for b in bonds_dict.values())
        sites_cations[c_lab] = {
            "element": c_el,
            "w": c_info["w"],
            "N": total_cn,
            "BVP": bvp_dict,
            "bonds": bonds_dict,
        }

        # Populate reciprocal N(Xi-M) in sites_anions using exact stoichiometric conservation:
        # w(M) * N(M-Xi) == w(Xi) * N(Xi-M)
        for x_lab, b_data in bonds_dict.items():
            w_X = sites_anions[x_lab]["w"]
            n_xm = int(round(c_info["w"] * b_data["N_MX"] / w_X))
            sites_anions[x_lab]["bonds"][c_lab] = n_xm

    for x_lab in sites_anions:
        sites_anions[x_lab]["N"] = sum(sites_anions[x_lab]["bonds"].values())

    # Find shortest cation-cation (acceptor) distances for MMCT and IVCT
    acceptors_list = []
    for c_lab in cat_labels:
        c_pos0 = expanded_orbits[c_lab]["coords"][0]
        for acc_lab in cat_labels:
            acc_el = expanded_orbits[acc_lab]["element"]
            acc_cn = sites_cations[acc_lab]["N"]
            min_d = 999.0
            for acc_pos in expanded_orbits[acc_lab]["coords"]:
                for sh in shifts:
                    d_vec = (acc_pos + sh) - c_pos0
                    dist = float(np.linalg.norm(M_cart @ d_vec))
                    if 1.5 < dist < min_d:
                        min_d = dist
            if min_d < 6.0 and get_acceptor_chi(acc_el, acc_cn) is not None:
                acceptors_list.append({
                    "doping_site": c_lab,
                    "acceptor_site": acc_lab,
                    "acceptor_element": acc_el,
                    "CN_prime": acc_cn,
                    "d_host": round(min_d, 4),
                })

    return CrystalStructureInput(
        name=name,
        formula=formula,
        space_group=space_group,
        Z=Z,
        V=V,
        sites_cations=sites_cations,
        sites_anions=sites_anions,
        acceptors=acceptors_list,
    )


def parse_cif_string(cif_text: str, filename: str = "Uploaded.cif") -> CrystalStructureInput:
    """
    Parses a crystallographic .cif file using pymatgen (if available) or the built-in pure-Python CIF parser.
    """
    # Extract cell parameters via regex (works on all standard CIFs)
    def _find_tag(tag: str, default: float = 0.0) -> float:
        m = re.search(rf"{tag}\s+([^\s#]+)", cif_text)
        return _parse_cif_float(m.group(1)) if m else default

    a = _find_tag(r"_cell_length_a", 5.0)
    b = _find_tag(r"_cell_length_b", a)
    c = _find_tag(r"_cell_length_c", a)
    alpha = _find_tag(r"_cell_angle_alpha", 90.0)
    beta = _find_tag(r"_cell_angle_beta", 90.0)
    gamma = _find_tag(r"_cell_angle_gamma", 90.0)
    Z = int(round(_find_tag(r"_cell_formula_units_Z", 1.0)))

    m_name = re.search(r"_chemical_formula_structural\s+([^\n\r]+)", cif_text) or re.search(r"_chemical_formula_sum\s+([^\n\r]+)", cif_text)
    formula = m_name.group(1).strip().strip("'").strip('"') if m_name else filename.replace(".cif", "")

    m_sg = re.search(r"_symmetry_space_group_name_H-M\s+([^\n\r]+)", cif_text) or re.search(r"_space_group_name_H-M_alt\s+([^\n\r]+)", cif_text)
    sg_num = int(round(_find_tag(r"_symmetry_Int_Tables_number", 0.0) or _find_tag(r"_space_group_IT_number", 0.0)))
    space_group = (m_sg.group(1).strip().strip("'").strip('"') + (f" ({sg_num})" if sg_num else "")) if m_sg else (str(sg_num) if sg_num else "P1")

    symops = []
    for m_op in re.finditer(r"['\"]([+-]?[xyz0-9/\.\s\+\-]+,[+-]?[xyz0-9/\.\s\+\-]+,[+-]?[xyz0-9/\.\s\+\-]+)['\"]", cif_text):
        symops.append(m_op.group(1).strip())

    # Parse _atom_site_ loop
    lines = [ln.strip() for ln in cif_text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    asymmetric_sites = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("loop_"):
            headers = []
            i += 1
            while i < len(lines) and lines[i].startswith("_"):
                headers.append(lines[i].split()[0])
                i += 1
            if any(h.startswith("_atom_site_fract_x") for h in headers):
                idx_x = headers.index("_atom_site_fract_x")
                idx_y = headers.index("_atom_site_fract_y")
                idx_z = headers.index("_atom_site_fract_z")
                idx_lbl = headers.index("_atom_site_label") if "_atom_site_label" in headers else None
                idx_sym = headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in headers else idx_lbl
                while i < len(lines) and not lines[i].startswith("loop_") and not lines[i].startswith("_"):
                    tokens = lines[i].split()
                    if len(tokens) >= len(headers):
                        lbl = tokens[idx_lbl] if idx_lbl is not None else tokens[idx_sym]
                        el = strip_element_symbol(tokens[idx_sym] if idx_sym is not None else lbl)
                        asymmetric_sites.append({
                            "label": lbl,
                            "element": el,
                            "x": tokens[idx_x],
                            "y": tokens[idx_y],
                            "z": tokens[idx_z],
                        })
                    i += 1
                continue
        i += 1

    return build_crystal_from_unit_cell_sites(
        name=filename,
        formula=formula,
        space_group=space_group,
        a=a, b=b, c=c,
        alpha_deg=alpha, beta_deg=beta, gamma_deg=gamma,
        Z=Z,
        asymmetric_sites=asymmetric_sites,
        symops=symops if symops else None,
    )


def inspect_json_for_missing_parameters(json_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Inspects an uploaded JSON structure (either Materials Project pymatgen format or custom
    summary format like mp-1227044.json) and identifies any missing lattice parameters so
    the GUI can prompt the user with intelligent pre-filled estimates.
    """
    # Case 1: Standard pymatgen structure JSON with 'lattice' and 'sites'
    if "lattice" in json_data and "sites" in json_data:
        lat = json_data["lattice"]
        return {
            "is_complete": True,
            "format": "pymatgen",
            "a": float(lat["a"]),
            "b": float(lat["b"]),
            "c": float(lat["c"]),
            "alpha": float(lat.get("alpha", 90.0)),
            "beta": float(lat.get("beta", 90.0)),
            "gamma": float(lat.get("gamma", 90.0)),
            "Z": int(json_data.get("Z", 2)),
            "space_group": json_data.get("space_group", "P1"),
            "missing_fields": [],
        }

    # Case 2: Custom summary JSON like mp-1227044.json ('crystal_structure' -> 'atomic_positions', 'symmetry', 'structure_meta')
    cs = json_data.get("crystal_structure", json_data)
    meta = cs.get("structure_meta", {})
    sym = cs.get("symmetry", {})
    positions = cs.get("atomic_positions", [])

    missing_fields = []
    lattice = cs.get("lattice", {})
    a = lattice.get("a")
    b = lattice.get("b")
    c = lattice.get("c")
    Z = meta.get("Z")

    if a is None:
        missing_fields.append("a (Lattice constant a in A)")
    if b is None:
        missing_fields.append("b (Lattice constant b in A)")
    if c is None:
        missing_fields.append("c (Lattice constant c in A)")
    if Z is None:
        missing_fields.append("Z (Formula units per unit cell)")

    # Estimate Z, Molar Mass, and Volume V from 'number_of_atoms' and 'density' if available!
    n_atoms_total = int(meta.get("number_of_atoms", 14))
    density_str = str(meta.get("density", "4.53"))
    m_rho = re.search(r"(\d+\.\d+|\d+)", density_str)
    rho = float(m_rho.group(1)) if m_rho else 4.53

    # Sum Wyckoff multiplicities from atomic_positions
    total_wyckoff = 0
    weighted_mass = 0.0
    formula_parts = []
    for pos in positions:
        w_str = str(pos.get("wyckoff", "2a"))
        m_w = re.match(r"^(\d+)", w_str)
        mult = int(m_w.group(1)) if m_w else 2
        el = strip_element_symbol(pos.get("element", "O"))
        total_wyckoff += mult
        weighted_mass += mult * ATOMIC_MASSES.get(el, 16.0)
        formula_parts.append((el, mult))

    # Deduce Z from GCD of Wyckoff multiplicities
    from math import gcd
    from functools import reduce
    mults = [m for _, m in formula_parts] if formula_parts else [2]
    est_Z = reduce(gcd, mults) if mults else 2
    formula_str = "".join(f"{el}{m//est_Z if (m//est_Z)>1 else ''}" for el, m in formula_parts)

    # Calculate exact unit cell volume V from density rho and total cell mass:
    # V (A^3) = (total_cell_mass_g_mol / (N_A * rho)) * 1e24 = total_cell_mass * 1.660539 / rho
    est_V = (weighted_mass * 1.6605390666) / rho if rho > 0 else 161.26

    # Estimate a, b, c based on crystal system
    cryst_sys = sym.get("crystal_system", "Tetragonal")
    if cryst_sys.lower() == "tetragonal":
        # For K2NiF4-type layered perovskites (c/a ~ 3.26 - 3.27)
        c_over_a = 3.273
        est_a = (est_V / c_over_a) ** (1.0 / 3.0)
        est_b = est_a
        est_c = est_a * c_over_a
    else:
        est_a = est_V ** (1.0 / 3.0)
        est_b = est_a
        est_c = est_a

    return {
        "is_complete": len(missing_fields) == 0,
        "format": "crystal_structure_json",
        "missing_fields": missing_fields,
        "estimated_Z": int(Z) if Z is not None else int(est_Z),
        "estimated_V": round(est_V, 2),
        "estimated_a": round(float(a) if a is not None else est_a, 4),
        "estimated_b": round(float(b) if b is not None else est_b, 4),
        "estimated_c": round(float(c) if c is not None else est_c, 4),
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "density_g_cm3": rho,
        "number_of_atoms": n_atoms_total,
        "space_group": f"{sym.get('symbol', 'I4mm')} (No. {sym.get('international_number', 107)})",
        "formula": formula_str or "CaYAlO4",
    }


def build_crystal_from_custom_json(
    json_data: Dict[str, Any],
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    Z: int,
    name: str = "Uploaded_JSON",
) -> CrystalStructureInput:
    """Builds a CrystalStructureInput from the uploaded JSON plus user-confirmed lattice parameters."""
    if "lattice" in json_data and "sites" in json_data:
        asym_sites = []
        for idx, s in enumerate(json_data["sites"]):
            el = s["species"][0]["element"]
            abc = s["abc"]
            asym_sites.append({"label": f"{el}{idx+1}", "element": el, "x": abc[0], "y": abc[1], "z": abc[2]})
        return build_crystal_from_unit_cell_sites(
            name=name, formula=json_data.get("formula", name),
            space_group=json_data.get("space_group", "P1"),
            a=a, b=b, c=c, alpha_deg=alpha, beta_deg=beta, gamma_deg=gamma, Z=Z,
            asymmetric_sites=asym_sites,
        )

    cs = json_data.get("crystal_structure", json_data)
    sym = cs.get("symmetry", {})
    sg_symbol = sym.get("symbol", "I4mm")
    positions = cs.get("atomic_positions", [])
    asym_sites = []
    el_counts: Dict[str, int] = {}
    for pos in positions:
        el = strip_element_symbol(pos["element"])
        el_counts[el] = el_counts.get(el, 0) + 1
        lbl = f"{el}{el_counts[el]}" if sum(1 for p in positions if strip_element_symbol(p["element"]) == el) > 1 else el
        asym_sites.append({
            "label": lbl,
            "element": el,
            "x": pos["a"],
            "y": pos["b"],
            "z": pos["c"],
        })

    info = inspect_json_for_missing_parameters(json_data)
    return build_crystal_from_unit_cell_sites(
        name=name,
        formula=info["formula"],
        space_group=sg_symbol,
        a=a, b=b, c=c,
        alpha_deg=alpha, beta_deg=beta, gamma_deg=gamma,
        Z=Z,
        asymmetric_sites=asym_sites,
    )


# ==============================================================================
# 4. BUILT-IN BENCHMARK DATASET FROM BOUTINAUD (2023) + CaYAlO4
# ==============================================================================

def get_boutinaud_2023_benchmarks() -> Dict[str, CrystalStructureInput]:
    """
    Returns 13 pre-configured crystal structures: 11 compounds (13 sites) from Tables 3, 4, and 5 of
    Boutinaud (PCCP 2023) plus CaYAlO4 in two structure models (ICDD 01-081-0742 and mp-1227044).
    """
    benchmarks: Dict[str, CrystalStructureInput] = {}

    # 1. CaYAlO4 (ICDD 01-081-0742 / Shannon 1992 Single-Crystal I4/mmm)
    benchmarks["CaYAlO4 (ICDD 01-081-0742, I4/mmm)"] = CrystalStructureInput(
        name="CaYAlO4 (ICDD 01-081-0742, Shannon 1992)",
        formula="CaYAlO4",
        space_group="I4/mmm (No. 139)",
        Z=2,
        V=157.77,
        sites_cations={
            "Y":  {"element": "Y",  "w": 2, "N": 9, "BVP": {"O1": 2.014, "O2": 2.014, "O3": 2.014},
                   "bonds": {"O2": {"N_MX": 1, "d_list": [2.247]}, "O3": {"N_MX": 4, "d_list": [2.459]*4}, "O1": {"N_MX": 4, "d_list": [2.592]*4}}},
            "Ca": {"element": "Ca", "w": 2, "N": 9, "BVP": {"O1": 1.967, "O2": 1.967, "O3": 1.967},
                   "bonds": {"O1": {"N_MX": 1, "d_list": [2.247]}, "O3": {"N_MX": 4, "d_list": [2.459]*4}, "O2": {"N_MX": 4, "d_list": [2.592]*4}}},
            "Al": {"element": "Al", "w": 2, "N": 6, "BVP": {"O1": 1.651, "O2": 1.651, "O3": 1.651},
                   "bonds": {"O3": {"N_MX": 4, "d_list": [1.823]*4}, "O2": {"N_MX": 1, "d_list": [1.995]}, "O1": {"N_MX": 1, "d_list": [1.995]}}},
        },
        sites_anions={
            "O1": {"element": "O", "w": 2, "N": 6, "bonds": {"Y": 4, "Ca": 1, "Al": 1}},
            "O2": {"element": "O", "w": 2, "N": 6, "bonds": {"Y": 1, "Ca": 4, "Al": 1}},
            "O3": {"element": "O", "w": 4, "N": 6, "bonds": {"Y": 2, "Ca": 2, "Al": 2}},
        },
        acceptors=[
            {"doping_site": "Y",  "acceptor_site": "Y", "acceptor_element": "Y",  "CN_prime": 9, "d_host": 3.626},
            {"doping_site": "Y",  "acceptor_site": "Bi_NNN", "acceptor_element": "Bi", "CN_prime": 9, "d_host": 5.155},
            {"doping_site": "Ca", "acceptor_site": "Y", "acceptor_element": "Y",  "CN_prime": 9, "d_host": 3.626},
            {"doping_site": "Ca", "acceptor_site": "Bi_NNN", "acceptor_element": "Bi", "CN_prime": 9, "d_host": 5.155},
        ],
        exp_EA_eV=4.43,  # 280 nm PLE peak
        exp_S_eV=1.31,   # 398 nm emission (3.12 eV) -> S = 1.31 eV
        reference_note="Experimental XRD (ICDD 01-081-0742). Excitation 280 nm, emission 398 / 420 nm (authors' unpublished data).",
    )

    # 2. CaYAlO4 (Materials Project mp-1227044 JSON, Tetragonal I4mm)
    benchmarks["CaYAlO4 (MP mp-1227044 JSON, I4mm)"] = CrystalStructureInput(
        name="CaYAlO4 (Materials Project mp-1227044, I4mm)",
        formula="CaYAlO4",
        space_group="I4mm (No. 107)",
        Z=2,
        V=161.26,
        sites_cations={
            "Y":  {"element": "Y",  "w": 2, "N": 9, "BVP": {"O1": 2.014, "O2": 2.014, "O3": 2.014},
                   "bonds": {"O2": {"N_MX": 1, "d_list": [2.2611]}, "O3": {"N_MX": 4, "d_list": [2.4550]*4}, "O1": {"N_MX": 4, "d_list": [2.6108]*4}}},
            "Ca": {"element": "Ca", "w": 2, "N": 9, "BVP": {"O1": 1.967, "O2": 1.967, "O3": 1.967},
                   "bonds": {"O1": {"N_MX": 1, "d_list": [2.2672]}, "O3": {"N_MX": 4, "d_list": [2.5611]*4}, "O2": {"N_MX": 4, "d_list": [2.6116]*4}}},
            "Al": {"element": "Al", "w": 2, "N": 6, "BVP": {"O1": 1.651, "O2": 1.651, "O3": 1.651},
                   "bonds": {"O3": {"N_MX": 4, "d_list": [1.8374]*4}, "O2": {"N_MX": 1, "d_list": [1.9783]}, "O1": {"N_MX": 1, "d_list": [2.0716]}}},
        },
        sites_anions={
            "O1": {"element": "O", "w": 2, "N": 6, "bonds": {"Y": 4, "Ca": 1, "Al": 1}},
            "O2": {"element": "O", "w": 2, "N": 6, "bonds": {"Y": 1, "Ca": 4, "Al": 1}},
            "O3": {"element": "O", "w": 4, "N": 6, "bonds": {"Y": 2, "Ca": 2, "Al": 2}},
        },
        acceptors=[
            {"doping_site": "Y",  "acceptor_site": "Y", "acceptor_element": "Y", "CN_prime": 9, "d_host": 3.666},
            {"doping_site": "Ca", "acceptor_site": "Y", "acceptor_element": "Y", "CN_prime": 9, "d_host": 3.626},
        ],
        exp_EA_eV=4.43,
        exp_S_eV=1.31,
        reference_note="DFT-relaxed I4mm structure from uploaded mp-1227044.json (density = 4.53 g/cm3, 14 atoms/cell).",
    )

    # 3. YPO4 (Zircon I41/amd, ICSD-79754, Boutinaud Table 5 & Table 6)
    benchmarks["YPO4 (ICSD-79754, Boutinaud Table 5 & 6)"] = CrystalStructureInput(
        name="YPO4 (Zircon I41/amd, ICSD-79754)",
        formula="YPO4",
        space_group="I41/amd (No. 141)",
        Z=4,
        V=284.98,
        sites_cations={
            "Y": {"element": "Y", "w": 4, "N": 8, "BVP": {"O": 2.014}, "bonds": {"O": {"N_MX": 8, "d_list": [2.313]*4 + [2.374]*4}}},
            "P": {"element": "P", "w": 4, "N": 4, "BVP": {"O": 1.617}, "bonds": {"O": {"N_MX": 4, "d_list": [1.535]*4}}},
        },
        sites_anions={"O": {"element": "O", "w": 16, "N": 3, "bonds": {"Y": 2, "P": 1}}},
        acceptors=[{"doping_site": "Y", "acceptor_site": "Y", "acceptor_element": "Y", "CN_prime": 8, "d_host": 3.766}],
        exp_EA_eV=5.48,
        exp_S_eV=0.38,
        pub_he_M={"Y": 1.12},
        pub_he_Bi={"Y": 1.22},
        reference_note="Published values: he(Y)=1.12, he(Bi3+)=1.22, EA=3.93 eV, EC=5.45 eV, ECT=4.95 eV. Stored excitation maximum 5.48 eV / S 0.38 eV: C band (5.46-5.51 eV, S 0.34-0.43 eV); the A band is not observed.",
    )

    # 4. Y3Al5O12 (YAG, Ia-3d, ICSD-280104, Boutinaud Table 5)
    benchmarks["Y3Al5O12 (YAG, ICSD-280104, Boutinaud Table 5)"] = CrystalStructureInput(
        name="Y3Al5O12 (Garnet Ia-3d, ICSD-280104)",
        formula="Y3Al5O12",
        space_group="Ia-3d (No. 230)",
        Z=8,
        V=1731.85,
        sites_cations={
            "Y":   {"element": "Y",  "w": 24, "N": 8, "BVP": {"O": 2.014}, "bonds": {"O": {"N_MX": 8, "d_list": [2.303]*4 + [2.438]*4}}},
            "Al1": {"element": "Al", "w": 16, "N": 6, "BVP": {"O": 1.651}, "bonds": {"O": {"N_MX": 6, "d_list": [1.921]*6}}},
            "Al2": {"element": "Al", "w": 24, "N": 4, "BVP": {"O": 1.651}, "bonds": {"O": {"N_MX": 4, "d_list": [1.761]*4}}},
        },
        sites_anions={"O": {"element": "O", "w": 96, "N": 4, "bonds": {"Y": 2, "Al1": 1, "Al2": 1}}},
        acceptors=[{"doping_site": "Y", "acceptor_site": "Y", "acceptor_element": "Y", "CN_prime": 8, "d_host": 3.677}],
        exp_EA_eV=4.55,
        exp_S_eV=0.53,
        pub_he_M={"Y": 0.94},
        pub_he_Bi={"Y": 1.02},
        reference_note="Published values: he(Y)=0.94, he(Bi3+)=1.02, E_exp=4.53-4.64 eV, S=0.49-0.58 eV, Delta_E=0.21 eV (A). Stored excitation maximum 4.55 eV / S 0.53 eV lie inside these ranges (A band).",
    )

    # 5. Y2O3 (Cubic Ia-3, ICSD-192862, S6 and C2 sites, Boutinaud Table 5)
    benchmarks["Y2O3 (ICSD-192862, S6 & C2 Sites, Boutinaud Table 5)"] = CrystalStructureInput(
        name="Y2O3 (Bixbyite Ia-3, ICSD-192862)",
        formula="Y2O3",
        space_group="Ia-3 (No. 206)",
        Z=16,
        V=1192.36,
        sites_cations={
            "Y1_S6": {"element": "Y", "w": 8,  "N": 6, "BVP": {"O": 2.014}, "bonds": {"O": {"N_MX": 6, "d_list": [2.284]*6}}},
            "Y2_C2": {"element": "Y", "w": 24, "N": 6, "BVP": {"O": 2.014}, "bonds": {"O": {"N_MX": 6, "d_list": [2.243, 2.243, 2.274, 2.274, 2.336, 2.336]}}},
        },
        sites_anions={"O": {"element": "O", "w": 48, "N": 4, "bonds": {"Y1_S6": 1, "Y2_C2": 3}}},
        acceptors=[
            {"doping_site": "Y1_S6", "acceptor_site": "Y2_C2", "acceptor_element": "Y", "CN_prime": 6, "d_host": 3.516},
            {"doping_site": "Y2_C2", "acceptor_site": "Y1_S6", "acceptor_element": "Y", "CN_prime": 6, "d_host": 3.516},
        ],
        exp_EA_eV=3.32,
        exp_S_eV=0.29,
        pub_he_M={"Y1_S6": 1.60, "Y2_C2": 1.61},
        pub_he_Bi={"Y1_S6": 1.72, "Y2_C2": 1.75},  # Table 5 (PDF p. 14): S6 1.60 / 1.72, C2 1.61 / 1.75
        reference_note="Published values: Y1[6]=S6 has he(M)=1.60, he(Bi)=1.72; Y2[6]=C2 has he(M)=1.61, he(Bi)=1.75.",
    )

    # 6. Cs2NaYCl6 (Elpasolite Fm-3m, ICSD-245353, Boutinaud Table 5)
    benchmarks["Cs2NaYCl6 (ICSD-245353, Boutinaud Table 5)"] = CrystalStructureInput(
        name="Cs2NaYCl6 (Elpasolite Fm-3m, ICSD-245353)",
        formula="Cs2NaYCl6",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=10.7315**3,
        sites_cations={
            "Y":  {"element": "Y",  "w": 4, "N": 6,  "BVP": {"Cl": 2.40}, "bonds": {"Cl": {"N_MX": 6,  "d_list": [2.6174]*6}}},
            "Na": {"element": "Na", "w": 4, "N": 6,  "BVP": {"Cl": 2.15}, "bonds": {"Cl": {"N_MX": 6,  "d_list": [2.7483]*6}}},
            "Cs": {"element": "Cs", "w": 8, "N": 12, "BVP": {"Cl": 2.79}, "bonds": {"Cl": {"N_MX": 12, "d_list": [3.7947]*12}}},
        },
        sites_anions={"Cl": {"element": "Cl", "w": 24, "N": 6, "bonds": {"Cs": 4, "Na": 1, "Y": 1}}},
        exp_EA_eV=3.84,
        exp_S_eV=0.12,
        pub_he_M={"Y": 1.24},
        pub_he_Bi={"Y": 1.32},
        reference_note="Published values: Cs2NaYCl6 (Y[6]): he(M)=1.24, he(Bi3+)=1.32, E_exp=3.84 eV.",
    )

    # 7. CaF2 (Fluorite Fm-3m, ICSD-656448, Boutinaud Table 4 & Table 5)
    benchmarks["CaF2 (ICSD-656448, Boutinaud Table 4 & 5)"] = CrystalStructureInput(
        name="CaF2 (Fluorite Fm-3m, ICSD-656448)",
        formula="CaF2",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=163.17,
        sites_cations={"Ca": {"element": "Ca", "w": 4, "N": 8, "BVP": {"F": 1.842}, "Q_override": 1.94, "bonds": {"F": {"N_MX": 8, "d_list": [2.3662]*8}}}},
        sites_anions={"F": {"element": "F", "w": 8, "N": 4, "bonds": {"Ca": 4}}},
        exp_EA_eV=5.86,
        exp_S_eV=0.35,
        pub_he_M={"Ca": 0.29},
        pub_he_Bi={"Ca": 0.29},
        reference_note="Published values: Q(Ca)=1.94, he(Ca)=0.29, he(Bi3+)=0.29, E_exp=5.86 eV.",
    )

    # 8. SrF2 (Fluorite Fm-3m, ICSD-404414, Boutinaud Table 5)
    benchmarks["SrF2 (ICSD-404414, Boutinaud Table 5)"] = CrystalStructureInput(
        name="SrF2 (Fluorite Fm-3m, ICSD-404414)",
        formula="SrF2",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=5.7996**3,
        sites_cations={"Sr": {"element": "Sr", "w": 4, "N": 8, "BVP": {"F": 2.019}, "bonds": {"F": {"N_MX": 8, "d_list": [2.5113]*8}}}},
        sites_anions={"F": {"element": "F", "w": 8, "N": 4, "bonds": {"Sr": 4}}},
        exp_EA_eV=5.78,
        exp_S_eV=0.31,
        pub_he_M={"Sr": 0.35},
        pub_he_Bi={"Sr": 0.32},
        reference_note="Published values: SrF2 (Sr[8]): he(M)=0.35, he(Bi3+)=0.32, E_exp=5.78 eV.",
    )

    # 9. BaF2 (Fluorite Fm-3m, ICSD-9009004, Boutinaud Table 5)
    benchmarks["BaF2 (ICSD-9009004, Boutinaud Table 5)"] = CrystalStructureInput(
        name="BaF2 (Fluorite Fm-3m, ICSD-9009004)",
        formula="BaF2",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=6.2001**3,
        sites_cations={"Ba": {"element": "Ba", "w": 4, "N": 8, "BVP": {"F": 2.188}, "bonds": {"F": {"N_MX": 8, "d_list": [2.6847]*8}}}},
        sites_anions={"F": {"element": "F", "w": 8, "N": 4, "bonds": {"Ba": 4}}},
        exp_EA_eV=5.75,
        exp_S_eV=1.00,
        pub_he_M={"Ba": 0.40},
        pub_he_Bi={"Ba": 0.33},
        reference_note="Published values: BaF2 (Ba[8]): he(M)=0.40, he(Bi3+)=0.33, E_exp=5.75 eV.",
    )

    # 10. CaO (Rock-salt Fm-3m, ICSD-7200686, Boutinaud Table 4 & Table 5)
    benchmarks["CaO (ICSD-7200686, Boutinaud Table 4 & 5)"] = CrystalStructureInput(
        name="CaO (Rock-salt Fm-3m, ICSD-7200686)",
        formula="CaO",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=111.33,
        sites_cations={"Ca": {"element": "Ca", "w": 4, "N": 6, "BVP": {"O": 1.967}, "Q_override": 1.835, "bonds": {"O": {"N_MX": 6, "d_list": [2.4053]*6}}}},
        sites_anions={"O": {"element": "O", "w": 4, "N": 6, "bonds": {"Ca": 6}}},
        exp_EA_eV=3.49,
        exp_S_eV=0.37,
        pub_he_M={"Ca": 1.35},
        pub_he_Bi={"Ca": 1.35},
        reference_note="Published values: Q(Ca)=1.835, he(Ca)=1.35, he(Bi3+)=1.35.",
    )

    # 11. CaS (Rock-salt Fm-3m, ICSD-603165, Boutinaud Table 4 & Table 5)
    benchmarks["CaS (ICSD-603165, Boutinaud Table 4 & 5)"] = CrystalStructureInput(
        name="CaS (Rock-salt Fm-3m, ICSD-603165)",
        formula="CaS",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=184.33,
        sites_cations={"Ca": {"element": "Ca", "w": 4, "N": 6, "BVP": {"S": 2.45}, "Q_override": 2.06, "bonds": {"S": {"N_MX": 6, "d_list": [2.8455]*6}}}},
        sites_anions={"S": {"element": "S", "w": 4, "N": 6, "bonds": {"Ca": 6}}},
        exp_EA_eV=3.01,
        exp_S_eV=0.27,
        pub_he_M={"Ca": 1.96},
        pub_he_Bi={"Ca": 1.96},
        reference_note="Published values: Q(Ca)=2.06, he(Ca)=1.96, he(Bi3+)=1.96.",
    )

    # 12. CaSe (Rock-salt Fm-3m, ICSD-41957, Boutinaud Table 4 & Table 5)
    benchmarks["CaSe (ICSD-41957, Boutinaud Table 4 & 5)"] = CrystalStructureInput(
        name="CaSe (Rock-salt Fm-3m, ICSD-41957)",
        formula="CaSe",
        space_group="Fm-3m (No. 225)",
        Z=4,
        V=206.92,
        sites_cations={"Ca": {"element": "Ca", "w": 4, "N": 6, "BVP": {"Se": 2.56}, "Q_override": 2.05, "bonds": {"Se": {"N_MX": 6, "d_list": [2.9574]*6}}}},
        sites_anions={"Se": {"element": "Se", "w": 4, "N": 6, "bonds": {"Ca": 6}}},
        exp_EA_eV=2.74,
        exp_S_eV=0.24,
        pub_he_M={"Ca": 2.10},
        pub_he_Bi={"Ca": 2.10},
        reference_note="Published values: Q(Ca)=2.05, he(Ca)=2.10, he(Bi3+)=2.10.",
    )

    # 13. CaTiO3 (Orthorhombic Pbnm, ICSD-165801, Boutinaud Table 3, 5 & 9)
    benchmarks["CaTiO3 (ICSD-165801, Boutinaud Table 3, 5 & 9)"] = CrystalStructureInput(
        name="CaTiO3 (Orthorhombic Perovskite Pbnm, ICSD-165801)",
        formula="CaTiO3",
        space_group="Pbnm (No. 62)",
        Z=4,
        V=223.62,
        sites_cations={
            "Ca": {"element": "Ca", "w": 4, "N": 8, "BVP": {"O1": 1.967, "O2": 1.967}, "Q_override": 1.93,
                   "bonds": {"O1": {"N_MX": 2, "d_list": [2.360, 2.463]}, "O2": {"N_MX": 6, "d_list": [2.382, 2.382, 2.625, 2.625, 2.671, 2.671]}}},
            "Ti": {"element": "Ti", "w": 4, "N": 6, "BVP": {"O1": 1.815, "O2": 1.815}, "Q_override": 4.11,
                   "bonds": {"O1": {"N_MX": 2, "d_list": [1.951, 1.951]}, "O2": {"N_MX": 4, "d_list": [1.956, 1.956, 1.961, 1.961]}}},
        },
        sites_anions={
            "O1": {"element": "O", "w": 4, "N": 4, "bonds": {"Ca": 2, "Ti": 2}},
            "O2": {"element": "O", "w": 8, "N": 5, "bonds": {"Ca": 3, "Ti": 2}},
        },
        acceptors=[{"doping_site": "Ca", "acceptor_site": "Ti", "acceptor_element": "Ti", "CN_prime": 6, "d_host": 3.169}],
        exp_EA_eV=3.76,
        exp_S_eV=1.55,
        pub_he_M={"Ca": 0.83, "Ti": 1.59},
        pub_he_Bi={"Ca": 0.83},
        reference_note="Published values: Ca[8], Ti[6] MMCT absorber & CT emitter (EA=4.41 eV, ECT=3.77 eV; Delta_E = 2.10-2.28 eV, CT type, and 2.28 eV in a second published table). Stored excitation maximum 3.76 eV / S 1.55 eV: CT-type band (3.65-3.87 eV, S 1.46-1.64 eV); the A-type value is 3.35 eV, S 1.13 eV.",
    )

    return benchmarks
