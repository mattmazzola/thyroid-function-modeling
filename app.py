# flake8: noqa
# pyright: reportMissingImports=false
import math
from typing import Optional

import streamlit as st  # type: ignore


# ---------------------------
# Helper functions
# ---------------------------

def convert_tsh_to_miu_per_l(value: float, unit: str) -> float:
    """Convert TSH to mIU/L (aka uIU/mL)."""
    unit = unit.strip().lower()
    if unit in {"miu/l", "mui/l", "uiu/ml", "µiu/ml", "uIU/mL".lower()}:
        return value
    if unit == "iu/l":
        return value * 1000.0
    # Fallback: assume mIU/L
    return value


def convert_concentration_to_mol_per_l(
    value: float,
    unit: str,
    molecular_weight_g_per_mol: float,
) -> float:
    """Convert concentration to mol/L using molecular weight when needed.

    Supports: pmol/L, nmol/L, pg/mL, ng/mL, ng/dL
    """
    unit = unit.strip().lower()
    if unit == "pmol/l":
        return value * 1e-12
    if unit == "nmol/l":
        return value * 1e-9
    if unit == "pg/ml":
        grams_per_l = value * 1e-12 / 1e-3  # pg -> g, mL -> L
        return grams_per_l / molecular_weight_g_per_mol
    if unit == "ng/ml":
        grams_per_l = value * 1e-9 / 1e-3
        return grams_per_l / molecular_weight_g_per_mol
    if unit == "ng/dl":
        grams_per_l = value * 1e-9 / 1e-1
        return grams_per_l / molecular_weight_g_per_mol
    # Fallback: assume pmol/L
    return value * 1e-12


def mol_per_l_to_pmol_per_l(value_mol_per_l: float) -> float:
    return value_mol_per_l * 1e12


def compute_spina_gt(
    tsh_miu_per_l: float,
    ft4_mol_per_l: float,
    *,
    alpha_t_l_inv: float = 0.1,
    beta_t_per_s: float = 1.1e-6,
    dt_miu_per_l: float = 2.75,
    k41_l_per_mol: float = 2.0e10,
    k42_l_per_mol: float = 2.0e8,
    tbg_mol_per_l: float = 300e-9,
    tbpa_mol_per_l: float = 4.5e-6,
) -> float:
    """Compute SPINA-GT (thyroid's secretory capacity) in mol/s.

    Formula (using free T4):
        GT = beta_T * (DT + [TSH]) / (alpha_T * [TSH]) * (1 + K41*[TBG] + K42*[TBPA]) * [FT4]

    Reference: Dietrich et al., Frontiers in Endocrinology 2016 (PMCID: PMC4899439)
    """
    if tsh_miu_per_l <= 0 or ft4_mol_per_l <= 0:
        return float("nan")
    protein_binding_term = 1.0 + k41_l_per_mol * tbg_mol_per_l + k42_l_per_mol * tbpa_mol_per_l
    gt_mol_per_s = (
        beta_t_per_s
        * (dt_miu_per_l + tsh_miu_per_l)
        / (alpha_t_l_inv * tsh_miu_per_l)
        * protein_binding_term
        * ft4_mol_per_l
    )
    return gt_mol_per_s


def compute_spina_gd(
    ft3_mol_per_l: float,
    ft4_mol_per_l: float,
    *,
    alpha31_l_inv: float = 0.026,
    beta31_per_s: float = 8.0e-6,
    km1_mol_per_l: float = 500e-9,  # 500 nmol/L
    k30_l_per_mol: float = 2.0e9,
    tbg_mol_per_l: float = 300e-9,
) -> float:
    """Compute SPINA-GD (sum activity of peripheral deiodinases) in mol/s.

    Formula (using free T3 and free T4):
        GD = beta31 * (KM1 + [FT4])/(1 + K30*[TBG]) * [FT3] / (alpha31 * [FT4])

    Reference: Dietrich et al., Frontiers in Endocrinology 2016 (PMCID: PMC4899439)
    """
    if ft3_mol_per_l <= 0 or ft4_mol_per_l <= 0:
        return float("nan")
    # For FT3-based calculation, reconstruct total T3 by multiplying with binding term
    # per SPINA references; omit this factor if using total T3 instead of FT3.
    binding_term = 1.0 + k30_l_per_mol * tbg_mol_per_l
    gd_mol_per_s = (
        beta31_per_s
        * (km1_mol_per_l + ft4_mol_per_l)
        * binding_term
        * (ft3_mol_per_l / (alpha31_l_inv * ft4_mol_per_l))
    )
    return gd_mol_per_s


def compute_tshi_jostel(tsh_miu_per_l: float, ft4_pmol_per_l: float, *, beta: float = 0.1345) -> float:
    """Compute Jostel's TSH Index (TSHI).

    TSHI = ln(TSH[mIU/L]) + 0.1345 * FT4[pmol/L]
    """
    if tsh_miu_per_l <= 0 or ft4_pmol_per_l <= 0:
        return float("nan")
    return math.log(tsh_miu_per_l) + beta * ft4_pmol_per_l


def ratio_ft3_rt3(ft3_mol_per_l: float, rt3_mol_per_l: float) -> Optional[float]:
    if ft3_mol_per_l > 0 and rt3_mol_per_l > 0:
        return ft3_mol_per_l / rt3_mol_per_l
    return None


# ---------------------------
# Streamlit UI
# ---------------------------

st.set_page_config(page_title="Thyroid Function Modeling", page_icon="🦋", layout="wide")
st.title("Thyroid Function Modeling")
st.caption("Enter lab values and adjust constants (optional) to calculate SPINA-GD, SPINA-GT, and Jostel's TSH index.")


# Custom styles for highlighted metric headings
st.markdown(
    """
    <style>
      .metric-heading { font-size: 1.4rem; font-weight: 700; margin: 0 0 0.25rem 0; }
      .metric-heading.gd { color: #1f77b4; }    /* blue */
      .metric-heading.gt { color: #2ca02c; }    /* green */
      .metric-heading.tshi { color: #d62728; }  /* red */
    </style>
    """,
    unsafe_allow_html=True,
)


with st.sidebar:
    st.header("Molecular weights and constants")
    st.caption("Defaults from Dietrich et al., Front Endocrinol 2016 and SPINA Thyr references.")

    # Molecular weights (g/mol)
    mw_t3 = st.number_input("Molecular weight T3 (g/mol)", value=650.97, min_value=1.0, step=0.01)
    mw_t4 = st.number_input("Molecular weight T4 (g/mol)", value=776.87, min_value=1.0, step=0.01)

    st.markdown("---")
    st.subheader("SPINA-GT constants")
    alpha_t = st.number_input("αT (L⁻¹)", value=0.1, min_value=0.0, step=0.001, format="%f")
    beta_t = st.number_input("βT (s⁻¹)", value=1.1e-6, min_value=0.0, step=1e-7, format="%e")
    dt = st.number_input("DT (mIU/L)", value=2.75, min_value=0.0, step=0.05)
    k41 = st.number_input("K41 (L/mol) T4–TBG", value=2.0e10, min_value=0.0, step=1e9, format="%e")
    k42 = st.number_input("K42 (L/mol) T4–TBPA", value=2.0e8, min_value=0.0, step=1e7, format="%e")
    tbg = st.number_input("[TBG] (mol/L)", value=300e-9, min_value=0.0, step=1e-9, format="%e")
    tbpa = st.number_input("[TBPA] (mol/L)", value=4.5e-6, min_value=0.0, step=1e-7, format="%e")

    st.markdown("---")
    st.subheader("SPINA-GD constants")
    alpha31 = st.number_input("α31 (L⁻¹)", value=0.026, min_value=0.0, step=0.001, format="%f")
    beta31 = st.number_input("β31 (s⁻¹)", value=8.0e-6, min_value=0.0, step=1e-6, format="%e")
    km1 = st.number_input("KM1 (mol/L) (500 nmol/L)", value=500e-9, min_value=0.0, step=1e-9, format="%e")
    k30 = st.number_input("K30 (L/mol) T3–TBG", value=2.0e9, min_value=0.0, step=1e8, format="%e")

    st.markdown("---")
    st.subheader("Jostel's TSH Index constant")
    beta_jti = st.number_input("β (dimensionless)", value=0.1345, min_value=0.0, step=0.0001, format="%f")


st.subheader("Lab inputs")

# Default units per spec
ft3_units = ["pg/mL", "pmol/L", "ng/dL"]
ft4_units = ["ng/dL", "pmol/L"]
rt3_units = ["ng/dL", "ng/mL", "pmol/L"]
tsh_units = ["uIU/mL", "mIU/L", "IU/L"]


def init_from_query_once():
    if st.session_state.get("_initialized_from_qp"):
        return
    # Pull from query params or use defaults
    qp = st.query_params
    st.session_state["ft3_val"] = float(qp.get("ft3", 3.2)) if "ft3" in qp else 3.2
    st.session_state["ft3_unit"] = str(qp.get("ft3_unit", "pg/mL"))

    st.session_state["ft4_val"] = float(qp.get("ft4", 1.2)) if "ft4" in qp else 1.2
    st.session_state["ft4_unit"] = str(qp.get("ft4_unit", "ng/dL"))

    st.session_state["rt3_val"] = float(qp.get("rt3", 14.0)) if "rt3" in qp else 14.0
    st.session_state["rt3_unit"] = str(qp.get("rt3_unit", "ng/dL"))

    st.session_state["tsh_val"] = float(qp.get("tsh", 2.0)) if "tsh" in qp else 2.0
    st.session_state["tsh_unit"] = str(qp.get("tsh_unit", "uIU/mL"))

    st.session_state["_initialized_from_qp"] = True


def update_url_from_state():
    st.query_params.update(
        {
            "ft3": f"{st.session_state['ft3_val']}",
            "ft3_unit": st.session_state["ft3_unit"],
            "ft4": f"{st.session_state['ft4_val']}",
            "ft4_unit": st.session_state["ft4_unit"],
            "rt3": f"{st.session_state['rt3_val']}",
            "rt3_unit": st.session_state["rt3_unit"],
            "tsh": f"{st.session_state['tsh_val']}",
            "tsh_unit": st.session_state["tsh_unit"],
        }
    )


# Initialize session_state once from the URL
init_from_query_once()


col1, col2, col3, col4 = st.columns(4)

with col1:
    st.number_input(
        "Free T3",
        min_value=0.0,
        step=0.1,
        key="ft3_val",
        on_change=update_url_from_state,
    )
    st.selectbox(
        "FT3 unit",
        options=ft3_units,
        key="ft3_unit",
        on_change=update_url_from_state,
    )

with col2:
    st.number_input(
        "Free T4",
        min_value=0.0,
        step=0.1,
        key="ft4_val",
        on_change=update_url_from_state,
    )
    st.selectbox(
        "FT4 unit",
        options=ft4_units,
        key="ft4_unit",
        on_change=update_url_from_state,
    )

with col3:
    st.number_input(
        "Reverse T3, serum",
        min_value=0.0,
        step=0.5,
        key="rt3_val",
        on_change=update_url_from_state,
    )
    st.selectbox(
        "rT3 unit",
        options=rt3_units,
        key="rt3_unit",
        on_change=update_url_from_state,
    )

with col4:
    st.number_input(
        "TSH",
        min_value=0.0,
        step=0.1,
        key="tsh_val",
        on_change=update_url_from_state,
    )
    st.selectbox(
        "TSH unit",
        options=tsh_units,
        key="tsh_unit",
        on_change=update_url_from_state,
    )


# Convert to SI for computations
ft3_mol_per_l = convert_concentration_to_mol_per_l(st.session_state["ft3_val"], st.session_state["ft3_unit"], mw_t3)
ft4_mol_per_l = convert_concentration_to_mol_per_l(st.session_state["ft4_val"], st.session_state["ft4_unit"], mw_t4)
rt3_mol_per_l = convert_concentration_to_mol_per_l(st.session_state["rt3_val"], st.session_state["rt3_unit"], mw_t3)
tsh_miu_per_l = convert_tsh_to_miu_per_l(st.session_state["tsh_val"], st.session_state["tsh_unit"])

# Convenience values in pmol/L for display and TSHI
ft4_pmol_per_l = mol_per_l_to_pmol_per_l(ft4_mol_per_l)
ft3_pmol_per_l = mol_per_l_to_pmol_per_l(ft3_mol_per_l)
rt3_pmol_per_l = mol_per_l_to_pmol_per_l(rt3_mol_per_l)


# Transparent calculations: unit conversions and formulas
with st.expander("Show how calculations are performed (unit conversions and formulas)"):
    st.markdown("### Unit conversions to pmol/L or mol/L used in formulas")
    # FT3 conversion explanation
    if st.session_state["ft3_unit"].lower() == "pg/ml":
        st.latex(
            rf"\mathrm{{FT3}}\,[pmol/L] = \mathrm{{FT3}}\,[pg/mL] \times \frac{{1000}}{{MW_{{T3}}}} = {st.session_state['ft3_val']:.3f} \times \frac{{1000}}{{{mw_t3:.2f}}} = {ft3_pmol_per_l:.2f}"
        )
    elif st.session_state["ft3_unit"].lower() == "ng/dl":
        st.latex(
            rf"\mathrm{{FT3}}\,[pmol/L] = \mathrm{{FT3}}\,[ng/dL] \times \frac{{10^4}}{{MW_{{T3}}}} = {st.session_state['ft3_val']:.3f} \times \frac{{10^4}}{{{mw_t3:.2f}}} = {ft3_pmol_per_l:.2f}"
        )
    else:
        st.latex(rf"\mathrm{{FT3}}\,[pmol/L] = {ft3_pmol_per_l:.2f}")

    # FT4 conversion explanation
    if st.session_state["ft4_unit"].lower() == "ng/dl":
        st.latex(
            rf"\mathrm{{FT4}}\,[pmol/L] = \mathrm{{FT4}}\,[ng/dL] \times \frac{{10^4}}{{MW_{{T4}}}} = {st.session_state['ft4_val']:.3f} \times \frac{{10^4}}{{{mw_t4:.2f}}} = {ft4_pmol_per_l:.2f}"
        )
    elif st.session_state["ft4_unit"].lower() == "pg/ml":
        st.latex(
            rf"\mathrm{{FT4}}\,[pmol/L] = \mathrm{{FT4}}\,[pg/mL] \times \frac{{1000}}{{MW_{{T4}}}} = {st.session_state['ft4_val']:.3f} \times \frac{{1000}}{{{mw_t4:.2f}}} = {ft4_pmol_per_l:.2f}"
        )
    else:
        st.latex(rf"\mathrm{{FT4}}\,[pmol/L] = {ft4_pmol_per_l:.2f}")

    # TSH note
    st.markdown(
        "- TSH units: 1 uIU/mL = 1 mIU/L; both are the same numerical value."
    )

    st.markdown("### SPINA expressions with live values")
    # Display SPINA-GT expression
    st.latex(
        r"\hat{G}_T=\beta_T\cdot \frac{D_T+ [TSH]}{\alpha_T\,[TSH]}\cdot (1+K_{41}[TBG]+K_{42}[TBPA])\cdot [FT4]"
    )
    st.latex(
        rf"= {beta_t:.2e} \cdot \frac{{{dt:.2f}+ {tsh_miu_per_l:.2f}}}{{{alpha_t:.3f}\times {tsh_miu_per_l:.2f}}}\cdot (1+{k41:.2e}\times {tbg:.2e}+{k42:.2e}\times {tbpa:.2e})\cdot {ft4_mol_per_l:.2e}\,\mathrm{{mol/s}}"
    )

    # Display SPINA-GD expression (FT3-based)
    st.latex(
        r"\hat{G}_D=\beta_{31}\cdot (K_M + [FT4])\cdot (1+K_{30}[TBG])\cdot \frac{[FT3]}{\alpha_{31}\,[FT4]}"
    )
    st.latex(
        rf"= {beta31:.2e}\cdot ({km1:.2e}+{ft4_mol_per_l:.2e})\cdot (1+{k30:.2e}\times {tbg:.2e})\cdot \frac{{{ft3_mol_per_l:.2e}}}{{{alpha31:.3f}\times {ft4_mol_per_l:.2e}}}\,\mathrm{{mol/s}}"
    )

    # Display TSHI expression (computed below in Results section)
    st.latex(r"TSHI = \ln([TSH]_{mIU/L}) + 0.1345\cdot [FT4]_{pmol/L}")
    st.latex(rf"= \ln({tsh_miu_per_l:.2f}) + {beta_jti:.4f}\times {ft4_pmol_per_l:.2f}")


st.markdown("---")
st.subheader("Results")

# Calculations
gt_mol_per_s = compute_spina_gt(
    tsh_miu_per_l,
    ft4_mol_per_l,
    alpha_t_l_inv=alpha_t,
    beta_t_per_s=beta_t,
    dt_miu_per_l=dt,
    k41_l_per_mol=k41,
    k42_l_per_mol=k42,
    tbg_mol_per_l=tbg,
    tbpa_mol_per_l=tbpa,
)
gd_mol_per_s = compute_spina_gd(
    ft3_mol_per_l,
    ft4_mol_per_l,
    alpha31_l_inv=alpha31,
    beta31_per_s=beta31,
    km1_mol_per_l=km1,
    k30_l_per_mol=k30,
    tbg_mol_per_l=tbg,
)
tshi = compute_tshi_jostel(tsh_miu_per_l, ft4_pmol_per_l, beta=beta_jti)

# Unit conversions for display
gt_pmol_per_s = gt_mol_per_s * 1e12 if math.isfinite(gt_mol_per_s) else float("nan")
gd_nmol_per_s = gd_mol_per_s * 1e9 if math.isfinite(gd_mol_per_s) else float("nan")


rcol1, rcol2, rcol3 = st.columns(3)

with rcol1:
    st.markdown('<div class="metric-heading gd">SPINA-GD</div>', unsafe_allow_html=True)
    st.metric(label="", value=f"{gd_nmol_per_s:,.2f} nmol/s" if math.isfinite(gd_nmol_per_s) else "—")
    with st.expander("What is SPINA-GD?"):
        st.write(
            """
            SPINA-GD (sum activity of peripheral deiodinases) estimates the maximum stimulated
            activity of step-up deiodination (conversion of T4 to T3). Typical reference range:
            ~20–60 nmol/s. Lower values suggest reduced conversion (e.g., NTIS, central hypothyroidism,
            certain medications); higher values may indicate hyperdeiodination.
            """
        )

with rcol2:
    st.markdown('<div class="metric-heading gt">SPINA-GT</div>', unsafe_allow_html=True)
    st.metric(label="", value=f"{gt_pmol_per_s:,.2f} pmol/s" if math.isfinite(gt_pmol_per_s) else "—")
    with st.expander("What is SPINA-GT?"):
        st.write(
            """
            SPINA-GT (thyroid's secretory capacity) estimates the thyroid’s maximum stimulated T4
            output. Adult reference range (70 kg calibration): 1.41–8.67 pmol/s. Low values point to
            reduced functional thyroid capacity (e.g., primary hypothyroidism, autoimmune thyroiditis);
            high values are seen in primary hyperthyroidism and goiter.
            """
        )

with rcol3:
    st.markdown('<div class="metric-heading tshi">Jostel\'s TSH Index</div>', unsafe_allow_html=True)
    st.metric(label="", value=f"{tshi:,.2f}" if math.isfinite(tshi) else "—")
    with st.expander("What is TSHI?"):
        st.write(
            """
            Jostel's TSH index (TSHI) is a marker of pituitary thyrotropic function based on the
            logarithmic standard model: TSHI = ln(TSH[mIU/L]) + 0.1345 × FT4[pmol/L]. Reference ~1.3–4.1.
            Lower values suggest thyrotropic insufficiency; higher values may reflect a higher central set point.
            """
        )


st.markdown("---")

st.subheader("Interpretation guidance")

ig1, ig2, ig3 = st.columns(3)
with ig1:
    st.markdown("**SPINA-GD (nmol/s)**")
    st.write(
        """
        - Low (< 20): reduced deiodination (e.g., NTIS, central hypothyroidism), certain drugs (e.g., propranolol).
        - Normal (20–60): typical euthyroid conversion capacity.
        - High (> 60): hyperdeiodination; may be seen in iodine deficiency or hyperthyroidism.
        """
    )
with ig2:
    st.markdown("**SPINA-GT (pmol/s)**")
    st.write(
        """
        - Low (< 1.41): reduced thyroid secretory capacity (primary hypothyroidism, untreated autoimmune thyroiditis).
        - Normal (1.41–8.67): typical adult range (calibrated for ~70 kg, plasma volume ~2.5 L).
        - High (> 8.67): primary hyperthyroidism (e.g., Graves’, toxic adenoma), diffuse/nodular goiter.
        """
    )
with ig3:
    st.markdown("**TSHI (unitless)**")
    st.write(
        """
        - Low (< 1.3): suggests thyrotropic insufficiency (secondary hypothyroidism); validate with clinical context.
        - Normal (1.3–4.1): typical central set point range.
        - High (> 4.1): higher central set point; may accompany stress/type 2 allostasis or subclinical hypothyroidism.
        """
    )


with st.expander("Additional ratios (optional)"):
    ft3_rt3 = ratio_ft3_rt3(ft3_mol_per_l, rt3_mol_per_l)
    colx, coly, colz = st.columns(3)
    with colx:
        st.write(f"FT3: {ft3_pmol_per_l:,.2f} pmol/L")
    with coly:
        st.write(f"FT4: {ft4_pmol_per_l:,.2f} pmol/L")
    with colz:
        st.write(f"rT3: {rt3_pmol_per_l:,.2f} pmol/L")
    if ft3_rt3 is not None:
        st.write(f"FT3/rT3 ratio: {ft3_rt3:,.2f}")
    else:
        st.write("FT3/rT3 ratio: —")


st.markdown("---")
st.caption(
    "This tool is for educational purposes and should not replace clinical judgment. References: Dietrich et al., Front Endocrinol (Lausanne) 2016; SPINA Thyr documentation."
)


