"""
bi3_format.py
=============
Typesetting helpers for the Bi3+ calculator (no Streamlit, no matplotlib): chemical formulas with
subscripts, Hermann-Mauguin space-group symbols, math-style table headers, and the live citation.

Every text function exists in up to three forms: *_html (web page), *_unicode (plain widgets such as
buttons and status messages) and *_math (matplotlib mathtext in the figures).
"""
import html
import json
import re
import time
import urllib.request

CITATION_URL = "https://omasgroup.org/citation2.json"
CITATION_TIMEOUT_S = 4
CITATION_REFRESH_S = 60  # a successful download is reused for this long
CITATION_RETRY_S = 20    # after a failure the program waits this long before asking the server again

_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_BAR = "̄"  # combining macron

# ------------------------------------------------------------------ chemical formulas
_NUM = r"\d+(?:\.\d+)?"


def _formula_tokens(formula):
    """Split a formula into (text, is_subscript) pieces.

    A number that follows an element symbol or a closing bracket is a subscript (Y2O3, Ca(OH)2);
    a number at the start or after a middle dot (CuSO4·5H2O) is a coefficient. An explicit count of 1
    (sum formulas such as 'Ca1 Y1 Al1 O4') is dropped and the blanks between elements are removed.
    """
    s = str(formula).strip().replace("*", "·")
    s = re.sub(r"\s+", "", s)
    out, i, prev = [], 0, ""  # prev: 'el', 'close', 'dot' or ''
    while i < len(s):
        ch = s[i]
        m = re.match(_NUM, s[i:])
        if m and ch.isdigit():
            num = m.group(0)
            if prev in ("el", "close"):
                if num not in ("1", "1.0"):
                    out.append((num, True))
            else:
                out.append((num, False))
            i += len(num)
            prev = "num"
            continue
        m = re.match(r"[A-Z][a-z]?", s[i:])
        if m:
            out.append((m.group(0), False))
            i += len(m.group(0))
            prev = "el"
            continue
        out.append((ch, False))
        prev = "close" if ch == ")" else ("dot" if ch == "·" else "")
        i += 1
    return out


def formula_html(formula):
    return "".join("<sub>%s</sub>" % html.escape(t) if sub else html.escape(t) for t, sub in _formula_tokens(formula))


def formula_unicode(formula):
    return "".join(t.translate(_SUB) if sub else t for t, sub in _formula_tokens(formula))


def formula_math(formula):
    """'Y3Al5O12' -> 'Y$_{3}$Al$_{5}$O$_{12}$' (real subscripts in matplotlib text)."""
    return "".join("$_{%s}$" % t if sub else t for t, sub in _formula_tokens(formula))


def label_parts(label):
    """Split 'Y2O3 (ICSD-192862, S6 & C2 Sites)' into ('Y2O3', ' (ICSD-192862, S6 & C2 Sites)')."""
    m = re.match(r"^(\S+)(.*)$", str(label).strip(), flags=re.S)
    return (m.group(1), m.group(2)) if m else (str(label), "")


def clean_benchmark_label(label):
    """Drop the literature table references from a built-in compound label."""
    s = re.sub(r",?\s*Boutinaud\s+Table\s+[\d,& ]+", "", str(label))
    return re.sub(r"\s+,", ",", s).replace("(,", "(").replace(" )", ")")


def compound_html(label):
    head, tail = label_parts(clean_benchmark_label(label))
    return formula_html(head) + html.escape(tail)


def compound_unicode(label):
    head, tail = label_parts(clean_benchmark_label(label))
    return formula_unicode(head) + tail


# ------------------------------------------------------------------ space groups
_SCREW = {"2": "1", "3": "12", "4": "123", "6": "12345"}
_SG_FIXED = {  # trigonal / cubic symbols whose digit pairs cannot be read from the string alone
    "P3112": ["P", "31", "1", "2"], "P3121": ["P", "31", "2", "1"],
    "P3212": ["P", "32", "1", "2"], "P3221": ["P", "32", "2", "1"],
    "P213": ["P", "21", "3"], "I213": ["I", "21", "3"], "P4232": ["P", "42", "3", "2"],
    "F4132": ["F", "41", "3", "2"], "I4132": ["I", "41", "3", "2"],
    "P4332": ["P", "43", "3", "2"], "P4132": ["P", "41", "3", "2"],
}
_SG_ITALIC = set("PIFABCRmabcdenghkl")


def split_space_group(sg):
    """'I4/mmm (No. 139)' -> ('I4/mmm', '139')."""
    s = str(sg).strip()
    m = re.match(r"^(.*?)\s*\((?:No\.?\s*)?(\d+)\)\s*$", s)
    if m:
        return m.group(1).strip(), m.group(2)
    return s, ""


def _sg_tokens(symbol):
    """Hermann-Mauguin symbol -> list of (kind, text): 'it' italic letter, 'n' number,
    'bar' number with an overbar, 'screw' (axis, subscript), 'txt' other characters."""
    s = symbol.strip().replace("_", "")  # '4_2' and '42' mean the same screw axis here
    if " " in symbol.strip() and "_" not in symbol:
        parts = symbol.strip().split()
    else:
        key = s.replace(" ", "")
        parts = _SG_FIXED.get(key)
        if parts is None:
            parts = []
            i = 0
            while i < len(key):
                c = key[i]
                if c == "-" and i + 1 < len(key) and key[i + 1].isdigit():
                    parts.append("-" + key[i + 1]); i += 2; continue
                if c.isdigit():
                    nxt = key[i + 1] if i + 1 < len(key) else ""
                    if c in _SCREW and nxt and nxt in _SCREW[c]:
                        parts.append(c + nxt); i += 2; continue
                parts.append(c); i += 1
    out = []
    for p in parts:
        if "/" in p and len(p) > 1:  # e.g. '4/mmm' given as one blank-free chunk
            p = p.replace("/", " / ")
            for q in p.split(" "):
                out.extend(_sg_piece(q) if q != "/" else [("txt", "/")])
        else:
            out.extend(_sg_piece(p))
    return out


def _sg_piece(p):
    out = []
    i = 0
    while i < len(p):
        c = p[i]
        if c == "-" and i + 1 < len(p) and p[i + 1].isdigit():
            out.append(("bar", p[i + 1])); i += 2
        elif c.isdigit():
            j = i
            while j < len(p) and p[j].isdigit():
                j += 1
            num = p[i:j]
            out.append(("screw", num) if len(num) == 2 else ("n", num))
            i = j
        elif c in _SG_ITALIC:
            out.append(("it", c)); i += 1
        else:
            out.append(("txt", c)); i += 1
    return out


def space_group_html(sg):
    symbol, number = split_space_group(sg)
    parts = []
    for kind, t in _sg_tokens(symbol):
        if kind == "it":
            parts.append("<i>%s</i>" % html.escape(t))
        elif kind == "bar":
            parts.append('<span style="text-decoration:overline">%s</span>' % t)
        elif kind == "screw":
            parts.append("%s<sub>%s</sub>" % (t[0], t[1]))
        else:
            parts.append(html.escape(t))
    s = "".join(parts).replace("</i><i>", "")
    return s + (" (No. %s)" % number if number else "")


def space_group_unicode(sg):
    symbol, number = split_space_group(sg)
    parts = []
    for kind, t in _sg_tokens(symbol):
        if kind == "bar":
            parts.append(t + _BAR)
        elif kind == "screw":
            parts.append(t[0] + t[1].translate(_SUB))
        else:
            parts.append(t)
    return "".join(parts) + (" (No. %s)" % number if number else "")


# ------------------------------------------------------------------ counts
def site_count_text(n):
    n = int(n)
    return "No sites" if n <= 0 else ("1 site" if n == 1 else "%d sites" % n)


# ------------------------------------------------------------------ math in table headers and cells
_HEADERS = {
    "Doping Site": "Doping site", "Raw Site": "Raw site", "CN": "CN", "Viable Bi3+ Site": "Viable Bi<sup>3+</sup> site",
    "Mean d_host (A)": "Mean <i>d</i><sub>host</sub> (Å)", "Delta_d (A)": "Δ<i>d</i> (Å)",
    "Q(M) [BVS]": "<i>Q</i>(M) (BVS)",
    "he(M) [Calc]": "<i>h</i><sub>e</sub>(M) calc.", "he(M) [Pub]": "<i>h</i><sub>e</sub>(M) publ.",
    "he(Bi3+) [Calc]": "<i>h</i><sub>e</sub>(Bi<sup>3+</sup>) calc.",
    "he(Bi3+) [Pub]": "<i>h</i><sub>e</sub>(Bi<sup>3+</sup>) publ.",
    "E_A (eV) [Eq. 18]": "<i>E</i><sub>A</sub> (eV)", "lambda_A (nm)": "<i>λ</i><sub>A</sub> (nm)",
    "E_A (cm^-1)": "<i>E</i><sub>A</sub> (cm<sup>−1</sup>)",
    "E_C (eV) [Eq. 19]": "<i>E</i><sub>C</sub> (eV)", "lambda_C (nm)": "<i>λ</i><sub>C</sub> (nm)",
    "Acceptor [CN']": "Acceptor [CN′]", "chi_CN'/d (A^-1)": "<i>χ</i><sub>CN′</sub>/<i>d</i> (Å<sup>−1</sup>)",
    "E_MMCT (eV) [Eq. 4]": "<i>E</i><sub>MMCT</sub> (eV)", "lambda_MMCT (nm)": "<i>λ</i><sub>MMCT</sub> (nm)",
    "E_IVCT Pair (eV) [Eq. 4]": "<i>E</i><sub>IVCT</sub> pair (eV)", "lambda_IVCT Pair (nm)": "<i>λ</i><sub>IVCT</sub> pair (nm)",
    "Exp. E_exc (eV)": "Exp. <i>E</i><sub>exc</sub> (eV)", "Exp. Stokes Shift S (eV)": "Exp. Stokes shift <i>S</i> (eV)",
    "Delta_E(Bi3+) (eV) [Eq. 5]": "Δ<i>E</i>(Bi<sup>3+</sup>) (eV)",
    "Excited-State Regime": "Excited-state regime", "Predicted Emission Origin": "Predicted emission origin",
    "CT parameters": "CT parameters",
    # Matrix I
    "Multiplicity w": "Multiplicity <i>w</i>", "Stoich. Coeff (w/Z)": "Stoich. coeff. (<i>w</i>/<i>Z</i>)",
    "Total CN (N)": "Total CN (<i>N</i>)", "BVS Q(M) [Eq. 10]": "BVS <i>Q</i>(M)",
    "Delta_d (A) [Table 1]": "Δ<i>d</i> (Å)",
    # sub-shell table
    "N(M-Xi)": "<i>N</i>(M–X<sub>i</sub>)", "N(Xi-M)": "<i>N</i>(X<sub>i</sub>–M)",
    "d_host (A)": "<i>d</i><sub>host</sub> (Å)", "d_Bi (A)": "<i>d</i><sub>Bi</sub> (Å)",
    "K(M-Xi) [Eq. 7]": "<i>K</i>(M–X<sub>i</sub>)", "K(Xi-M) [Eq. 7]": "<i>K</i>(X<sub>i</sub>–M)",
    "Q(Xi) [Eq. 9]": "<i>Q</i>(X<sub>i</sub>)", "N*(M-Xi) [Eq. 15]": "<i>N</i>*(M–X<sub>i</sub>)",
    "v_host (A^3) [Eq. 16]": "<i>v</i><sub>host</sub> (Å<sup>3</sup>)",
    "Gamma_1 (host) [Eq. 12]": "Γ<sub>1</sub> (host)", "Gamma_2 [Eq. 13]": "Γ<sub>2</sub>",
    "Gamma_3 (host) [Eq. 14]": "Γ<sub>3</sub> (host)", "f_c (host) [Eq. 11]": "<i>f</i><sub>c</sub> (host)",
    "alpha (host, A^3) [Eq. 17]": "<i>α</i> (host, Å<sup>3</sup>)",
    "f_c (Bi3+)": "<i>f</i><sub>c</sub> (Bi<sup>3+</sup>)", "alpha (Bi3+, A^3)": "<i>α</i> (Bi<sup>3+</sup>, Å<sup>3</sup>)",
    # benchmark table
    "Benchmark Compound": "Compound", "Site [CN]": "Site [CN]",
    "Calc he(M)": "<i>h</i><sub>e</sub>(M) calc.", "Pub he(M)": "<i>h</i><sub>e</sub>(M) publ.",
    "Calc he(Bi3+)": "<i>h</i><sub>e</sub>(Bi<sup>3+</sup>) calc.", "Pub he(Bi3+)": "<i>h</i><sub>e</sub>(Bi<sup>3+</sup>) publ.",
    "Calc E_A (eV)": "<i>E</i><sub>A</sub> calc. (eV)", "Calc E_C (eV)": "<i>E</i><sub>C</sub> calc. (eV)",
    "Calc E_MMCT (eV)": "<i>E</i><sub>MMCT</sub> calc. (eV)", "Exp E_exc (eV)": "Exp. <i>E</i><sub>exc</sub> (eV)",
    "Calc Delta_E (eV)": "Δ<i>E</i> calc. (eV)", "Emission Origin": "Emission origin",
    "Exp. band type": "Exp. band type", "Exp. value from": "Exp. value from",
}


def plain_column(name):
    """Column name for CSV files: the same name without the internal equation / table tags."""
    return re.sub(r"\s*\[(?:Eq\.[^\]]*|Table[^\]]*)\]", "", str(name)).strip()


def plain_columns(df):
    return df.rename(columns={c: plain_column(c) for c in df.columns})


def math_html(text):
    """Escape a text and typeset the symbols the program uses (E_A, Delta_E, Bi3+, 3P1 -> 1S0, Angstrom)."""
    s = html.escape(str(text), quote=False)
    s = re.sub(r"Delta_E", "Δ<i>E</i>", s)
    s = re.sub(r"\bE_(CT|MMCT|IVCT|exc|A|C)\b", r"<i>E</i><sub>\1</sub>", s)
    s = re.sub(r"\bd_host\b", "<i>d</i><sub>host</sub>", s)
    s = re.sub(r"\bd_Bi\b", "<i>d</i><sub>Bi</sub>", s)
    s = s.replace("Bi3+", "Bi<sup>3+</sup>")
    s = re.sub(r"\b(\d)([SPD])(\d)\b", r"<sup>\1</sup>\2<sub>\3</sub>", s)
    s = s.replace("-&gt;", "→").replace("&lt;=", "≤").replace("&gt;=", "≥")
    s = re.sub(r"(\d) A\)", r"\1 Å)", s)
    return s


def header_html(col):
    col = str(col)
    if col in _HEADERS:
        return _HEADERS[col]
    m = re.match(r"^N\(Site-(.+)\)$", col)
    if m:
        return "<i>N</i>(Site–%s)" % html.escape(m.group(1))
    return math_html(plain_column(col))


def _cell(v):
    try:
        import pandas as pd
        if v is None or (not isinstance(v, str) and pd.isna(v)):
            return ""
    except Exception:
        if v is None:
            return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return format(v, ".10g")
    return math_html(v)


def table_html(df, compound_col=None):
    """A DataFrame as an HTML table with typeset headers (scrolls sideways when wide)."""
    head = "".join("<th>%s</th>" % header_html(c) for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            num = isinstance(v, (int, float)) and not isinstance(v, bool)
            txt = compound_html(v) if c == compound_col else _cell(v)
            cells.append('<td class="%s">%s</td>' % ("n" if num else "t", txt))
        rows.append("<tr>%s</tr>" % "".join(cells))
    return ('<div class="bi3tab"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (head, "".join(rows)))


TABLE_CSS = """
<style>
.bi3tab{overflow-x:auto;margin:0 0 1rem 0;border:1px solid rgba(128,128,128,.35);border-radius:6px}
.bi3tab table{border-collapse:collapse;width:100%;font-size:.85rem}
.bi3tab th{background:rgba(128,128,128,.16);text-align:left;padding:6px 10px;white-space:nowrap;
 border-bottom:1px solid rgba(128,128,128,.5);font-weight:600}
.bi3tab td{padding:5px 10px;border-bottom:1px solid rgba(128,128,128,.18);white-space:nowrap;vertical-align:top}
.bi3tab td.n{text-align:right;font-variant-numeric:tabular-nums}
.bi3tab tr:last-child td{border-bottom:none}
.bi3card{border:1px solid rgba(128,128,128,.35);border-radius:8px;padding:10px 14px;height:100%}
.bi3card .l{font-size:.8rem;opacity:.7;margin-bottom:2px}
.bi3card .v{font-size:1.45rem;line-height:1.3}
</style>
"""


def card_html(label, value_html):
    return '<div class="bi3card"><div class="l">%s</div><div class="v">%s</div></div>' % (label, value_html)


# ------------------------------------------------------------------ live citation
_cache = {"t": 0.0, "data": None, "fail_t": 0.0, "err": ""}


def fetch_citation(url=CITATION_URL, force=False):
    """Download the citation record. Returns (dict or None, error text).

    The file is read from the web each time it is needed; a successful download is reused for
    CITATION_REFRESH_S seconds so that clicking in the program does not repeat the request, and after
    a failure the program does not ask again for CITATION_RETRY_S seconds (so that an offline PC is
    not slowed down). Nothing is stored on disk and a failed download is never replaced by an old copy.
    """
    now = time.time()
    if not force and _cache["data"] is not None and now - _cache["t"] < CITATION_REFRESH_S:
        return _cache["data"], ""
    if not force and _cache["data"] is None and now - _cache["fail_t"] < CITATION_RETRY_S:
        return None, _cache["err"]
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Bi3-Calculator", "Cache-Control": "no-cache"})
        with urllib.request.urlopen(req, timeout=CITATION_TIMEOUT_S) as resp:
            data = json.loads(resp.read().decode("utf-8-sig"))
        if not isinstance(data, dict):
            raise ValueError("the file is not a JSON object")
        _cache.update(t=now, data=data)
        return data, ""
    except Exception as e:  # no internet, blocked, bad JSON ...
        _cache.update(t=0.0, data=None, fail_t=now, err=str(e))
        return None, str(e)


def citation_text(data):
    """One reference line from the fields authors, title, journal, volume, pages, year, doi, note."""
    g = lambda k: str(data.get(k, "") or "").strip()
    head = ", ".join(x for x in (g("authors"), g("title")) if x)
    jour = " ".join(x for x in (g("journal"), g("volume")) if x)
    tail = ", ".join(x for x in (jour, g("pages")) if x)
    if g("year"):
        tail = (tail + " " if tail else "") + "(%s)" % g("year")
    text = ", ".join(x for x in (head, tail) if x)
    if g("doi"):
        text += ", https://doi.org/" + g("doi")
    if g("note"):
        text += ". " + g("note")
    return text


def citation_for_report():
    """Text for the 'How to cite' block of the report (live download; says so when it fails)."""
    data, err = fetch_citation()
    if data is None:
        return "Could not be loaded (no internet connection?). Please read the current reference at %s" % CITATION_URL
    return citation_text(data)
