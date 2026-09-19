import os
import sys
import shutil
import subprocess
import time
import math
import tempfile
import ast
import numpy as np
import pandas as pd
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse

BASE = os.path.dirname(os.path.abspath(__file__))
VALIDATION_FILE = os.path.join(BASE, "GenCon_Validation_361.csv")
PREDICTION_FILE = os.path.join(BASE, "GenCon_Validation_361_Predictions.csv")
AMRFINDER_FILE = os.path.join(BASE, "GenCon_361_AMRFinder_Combined.tsv")

KEY_GENES = ["blaOXA-23", "blaOXA-24", "blaOXA-237", "blaOXA-72", "blaKPC-2"]

def normalize_genome_id(value):
    if pd.isna(value): return ""
    val = str(value).strip()
    if val.endswith(".0"):
        try: val = str(int(float(val)))
        except: pass
    return val

def parse_gene_list(value):
    if value is None or pd.isna(value): return []
    text = str(value).strip()
    if not text or text.lower() == "nan": return []
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, list):
            return sorted(set([str(x).split("=")[0].strip() for x in parsed if str(x).split("=")[0].strip()]))
    except:
        pass
    genes = []
    for item in text.split(","):
        gene = item.strip().split("=")[0].strip()
        if gene: genes.append(gene)
    return sorted(set(genes))

def get_gencon_genes(row):
    if "GenCon_Genes_Detected" in row.index:
        g = parse_gene_list(row["GenCon_Genes_Detected"])
        if g: return g
    for col in ["Gene_list", "AMR genotypes"]:
        if col in row.index:
            g = parse_gene_list(row[col])
            sel = [x for x in g if x in KEY_GENES]
            if sel: return sorted(set(sel))
    sel = []
    for gene in KEY_GENES:
        if gene in row.index and str(row[gene]).strip().upper() in ["1", "TRUE", "YES", "PRESENT", "COMPLETE"]:
            sel.append(gene)
    return sorted(set(sel))

validation_df = pd.DataFrame()
amrfinder_df = pd.DataFrame()

if os.path.isfile(PREDICTION_FILE):
    validation_df = pd.read_csv(PREDICTION_FILE)
elif os.path.isfile(VALIDATION_FILE):
    validation_df = pd.read_csv(VALIDATION_FILE)

if not validation_df.empty:
    if "BV_BRC_genome_id" in validation_df.columns:
        validation_df["BV_BRC_genome_id"] = validation_df["BV_BRC_genome_id"].apply(normalize_genome_id)
    if "Meropenem_AST" not in validation_df.columns:
        validation_df["Meropenem_AST"] = validation_df.get("Meropenem", "")
    validation_df["Meropenem_AST"] = validation_df["Meropenem_AST"].fillna("").astype(str).str.strip().str.upper()
    validation_df["GenCon_Genes_Detected"] = validation_df.apply(get_gencon_genes, axis=1)
    validation_df["GenCon_Prediction"] = validation_df["GenCon_Genes_Detected"].apply(lambda x: "R" if len(x) > 0 else "S")
    validation_df["Status"] = validation_df.apply(
        lambda r: "AST NOT PROVIDED" if r["Meropenem_AST"] not in ["R", "S"]
        else ("CONCORDANT" if r["GenCon_Prediction"] == r["Meropenem_AST"] else "DISCORDANT"),
        axis=1
    )

if os.path.isfile(AMRFINDER_FILE):
    try:
        amrfinder_df = pd.read_csv(AMRFINDER_FILE, sep="\t")
        if "BV_BRC_genome_id" in amrfinder_df.columns:
            amrfinder_df["BV_BRC_genome_id"] = amrfinder_df["BV_BRC_genome_id"].apply(normalize_genome_id)
    except Exception:
        pass

def safe_div(a, b): return (a / b) if b != 0 else 0.0

metrics = {
    "total": int(len(validation_df)), "actual_R": 0, "actual_S": 0, "predicted_R": 0, "predicted_S": 0,
    "tp": 0, "tn": 0, "fp": 0, "fn": 0, "accuracy": 0.0, "sensitivity": 0.0,
    "specificity": 0.0, "precision": 0.0, "f1": 0.0, "mcc": 0.0, "concordant": 0, "discordant": 0
}

if not validation_df.empty:
    actual = validation_df["Meropenem_AST"]
    predicted = validation_df["GenCon_Prediction"]
    valid = actual.isin(["R", "S"])
    act_v, pred_v = actual[valid], predicted[valid]
    tp = int(((act_v == "R") & (pred_v == "R")).sum())
    tn = int(((act_v == "S") & (pred_v == "S")).sum())
    fp = int(((act_v == "S") & (pred_v == "R")).sum())
    fn = int(((act_v == "R") & (pred_v == "S")).sum())
    total = tp + tn + fp + fn
    prec = safe_div(tp, tp + fp)
    sens = safe_div(tp, tp + fn)
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    metrics.update({
        "actual_R": int((actual == "R").sum()),
        "actual_S": int((actual == "S").sum()),
        "predicted_R": int((predicted == "R").sum()),
        "predicted_S": int((predicted == "S").sum()),
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "accuracy": safe_div(tp + tn, total),
        "sensitivity": sens,
        "specificity": safe_div(tn, tn + fp),
        "precision": prec,
        "f1": safe_div(2 * prec * sens, prec + sens),
        "mcc": safe_div((tp * tn) - (fp * fn), denom),
        "concordant": int((validation_df["Status"] == "CONCORDANT").sum()),
        "discordant": int((validation_df["Status"] == "DISCORDANT").sum())
    })

def clean_json_value(v):
    if pd.isna(v): return None
    if isinstance(v, (np.integer, int)): return int(v)
    if isinstance(v, (np.floating, float)): return float(v)
    return v

def row_to_json(row):
    out = {str(k): clean_json_value(v) for k, v in row.items()}
    if "GenCon_Genes_Detected" in out:
        out["GenCon_Genes_Detected"] = parse_gene_list(out["GenCon_Genes_Detected"])
    return out

def run_fasta_analysis(fasta_path, ast_result=""):
    amrfinder = shutil.which("amrfinder") or "/usr/local/bin/amrfinder"
    if not os.path.exists(amrfinder):
        return {"success": False, "message": "AMRFinderPlus binary not found."}

    out_file = os.path.join(tempfile.gettempdir(), f"res_{os.getpid()}_{int(time.time())}.tsv")
    cmd = [amrfinder, "-n", fasta_path, "-O", "Acinetobacter_baumannii", "-o", out_file]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            return {"success": False, "message": proc.stderr or "AMRFinder analysis failed."}
        res_df = pd.read_csv(out_file, sep="\t")
    except Exception as e:
        return {"success": False, "message": str(e)}
    finally:
        if os.path.exists(out_file):
            os.remove(out_file)

    detected = sorted(set(res_df["Element symbol"].dropna().astype(str).str.strip())) if "Element symbol" in res_df.columns else []
    gencon_detected = sorted(set(detected).intersection(KEY_GENES))
    prediction = "R" if gencon_detected else "S"
    ast_val = str(ast_result or "").strip().upper()
    status = "AST NOT PROVIDED" if ast_val not in ["R", "S"] else ("CONCORDANT" if ast_val == prediction else "DISCORDANT")

    cols = [c for c in ["Element symbol", "Element name", "Class", "Method"] if c in res_df.columns]
    return {
        "success": True,
        "prediction": prediction,
        "prediction_text": "RESISTANT" if prediction == "R" else "SUSCEPTIBLE",
        "gencon_genes": gencon_detected,
        "all_detected_elements": detected,
        "ast": ast_val or "Not provided",
        "status": status,
        "reason": "At least one GenCon-AMR resistance determinant detected." if gencon_detected else "None of the GenCon-AMR resistance determinants detected.",
        "results": res_df[cols].replace({np.nan: None}).to_dict(orient="records")
    }

app = FastAPI(title="GenCon-AMR")

HTML_PAGE = """<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GenCon-AMR</title>
<style>
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f4f7fb; color: #172033; }
.navbar { background: #0d1b2a; color: white; padding: 15px 5%; display: flex; justify-content: space-between; align-items: center; }
.logo { font-size: 20px; font-weight: 800; } .logo span { color: #38bdf8; }
.nav-btn { background: transparent; border: 0; color: #dbeafe; padding: 8px 12px; cursor: pointer; font-size: 14px; }
.page { display: none; max-width: 1200px; margin: auto; padding: 30px 5%; } .page.active { display: block; }
.hero { background: linear-gradient(135deg, #0f2742, #123b59); color: white; border-radius: 18px; padding: 40px; margin-bottom: 25px; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin-bottom: 20px; }
.card { background: white; border-radius: 14px; padding: 20px; border: 1px solid #e3e9f1; }
.metric { font-size: 30px; font-weight: 800; } .metric-label { font-size: 12px; color: #667085; }
.btn { border: none; border-radius: 8px; padding: 10px 16px; cursor: pointer; font-weight: 700; }
.btn-primary { background: #0ea5e9; color: white; }
.btn-secondary { background: #e8f1f8; color: #123b59; }
.input { width: 100%; padding: 10px; border: 1px solid #d7dee8; border-radius: 8px; margin: 6px 0 12px; }
.upload-box { border: 2px dashed #9db8cb; border-radius: 10px; padding: 20px; background: #f8fbfd; margin: 10px 0; }
.prediction { border-radius: 12px; padding: 20px; text-align: center; margin-top: 15px; }
.prediction.r { background: #fff1f1; border: 1px solid #fecaca; } .prediction.s { background: #effdf4; border: 1px solid #bbf7d0; }
.gene { display: inline-block; background: #edf5fa; color: #0f4c66; border-radius: 6px; padding: 4px 8px; margin: 3px; font-family: monospace; }
table { width: 100%; border-collapse: collapse; margin-top: 15px; background: white; }
th, td { padding: 10px; border-top: 1px solid #e8edf3; text-align: left; font-size: 13px; }
th { background: #edf2f7; }
.r-text { color: #b91c1c; font-weight: 800; } .s-text { color: #166534; font-weight: 800; }
.rule { font-family: monospace; background: #0f172a; color: #e2e8f0; padding: 16px; border-radius: 10px; line-height: 1.7; overflow-x: auto; }
</style>
</head>
<body>
<div class="navbar">
    <div class="logo">🧬 GenCon<span>-AMR</span></div>
    <div>
        <button class="nav-btn" onclick="showPage('home')">Home</button>
        <button class="nav-btn" onclick="showPage('analysis')">Analyze Genome</button>
        <button class="nav-btn" onclick="showPage('dataset')">Validation Dataset</button>
        <button class="nav-btn" onclick="showPage('methods')">Methods</button>
    </div>
</div>

<section id="home" class="page active">
    <div class="hero">
        <h1>GenCon-AMR</h1>
        <p>Genotype–Phenotype Discordance Analyzer for Meropenem Resistance in <i>Acinetobacter baumannii</i>.</p>
        <button class="btn btn-primary" onclick="showPage('analysis')">🧬 Analyze a Genome</button>
    </div>
    <div id="homeMetrics" class="grid"></div>
</section>

<section id="analysis" class="page">
    <h2>🧬 Analyze Genome</h2>
    <div class="card">
        <label><b>Genome FASTA (.fasta, .fa, .fna)</b></label>
        <div class="upload-box"><input id="fastaFile" type="file" accept=".fasta,.fa,.fna"></div>
        <label><b>Known Meropenem AST (Optional)</b></label>
        <select id="fastaAST" class="input">
            <option value="">Not provided</option>
            <option value="R">Resistant (R)</option>
            <option value="S">Susceptible (S)</option>
        </select>
        <button class="btn btn-primary" onclick="analyzeGenome()">🔬 Run GenCon-AMR</button>
        <div id="loader" style="display:none; padding:12px; color:#0284c7; font-weight:bold;">Running AMRFinderPlus... Please wait.</div>
        <div id="analysisResult"></div>
    </div>
</section>

<section id="dataset" class="page">
    <h2>📊 Validation Dataset</h2>
    <div class="card">
        <input id="search" class="input" placeholder="Search Genome ID..." oninput="loadDataset()">
        <div id="datasetTable"></div>
    </div>
</section>

<section id="methods" class="page">
    <h2>ℹ️ Methods</h2>
    <div class="card">
        <div class="rule">
IF blaOXA-23 detected   → RESISTANT (R)<br>
OR blaOXA-24 detected   → RESISTANT (R)<br>
OR blaOXA-237 detected  → RESISTANT (R)<br>
OR blaOXA-72 detected   → RESISTANT (R)<br>
OR blaKPC-2 detected    → RESISTANT (R)<br>
OTHERWISE               → SUSCEPTIBLE (S)
        </div>
    </div>
</section>

<script>
function showPage(id) {
    document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
    document.getElementById(id).classList.add('active');
    if (id === 'home') loadSummary();
    if (id === 'dataset') loadDataset();
}

async function loadSummary() {
    const res = await fetch("/api/summary");
    const d = await res.json();
    document.getElementById("homeMetrics").innerHTML = `
        <div class="card"><div class="metric">${d.total}</div><div class="metric-label">Total Isolates</div></div>
        <div class="card"><div class="metric r-text">${d.actual_R}</div><div class="metric-label">AST Resistant</div></div>
        <div class="card"><div class="metric s-text">${d.actual_S}</div><div class="metric-label">AST Susceptible</div></div>
        <div class="card"><div class="metric">${d.accuracy}%</div><div class="metric-label">Accuracy</div></div>
        <div class="card"><div class="metric">${d.discordant}</div><div class="metric-label">Discordant</div></div>
    `;
}

async function loadDataset() {
    const q = encodeURIComponent(document.getElementById("search") ? document.getElementById("search").value : "");
    const res = await fetch(`/api/dataset?search=${q}`);
    const data = await res.json();
    let rows = data.map(r => `
        <tr>
            <td><b>${r.BV_BRC_genome_id}</b></td>
            <td><span class="${r.Meropenem_AST === 'R' ? 'r-text':'s-text'}">${r.Meropenem_AST || '—'}</span></td>
            <td><span class="${r.GenCon_Prediction === 'R' ? 'r-text':'s-text'}">${r.GenCon_Prediction}</span></td>
            <td>${r.Status}</td>
            <td>${(r.GenCon_Genes_Detected || []).join(", ") || "None"}</td>
        </tr>
    `).join('');
    document.getElementById("datasetTable").innerHTML = `
        <div style="overflow-x:auto;"><table><thead><tr><th>Genome ID</th><th>AST</th><th>GenCon</th><th>Status</th><th>Determinants</th></tr></thead><tbody>${rows}</tbody></table></div>
    `;
}

async function analyzeGenome() {
    const file = document.getElementById("fastaFile").files[0];
    const ast = document.getElementById("fastaAST").value;
    if (!file) return alert("Select a FASTA file first.");
    document.getElementById("loader").style.display = "block";
    document.getElementById("analysisResult").innerHTML = "";

    const fd = new FormData();
    fd.append("file", file);
    fd.append("ast", ast);

    try {
        const res = await fetch("/api/analyze-fasta", { method: "POST", body: fd });
        const d = await res.json();
        if (!d.success) {
            document.getElementById("analysisResult").innerHTML = `<div style="color:red; margin-top:10px;"><b>Error:</b> ${d.message}</div>`;
            return;
        }
        const cls = d.prediction === "R" ? "r" : "s";
        const genes = d.gencon_genes.map(g => `<span class="gene">${g}</span>`).join("") || "None";
        document.getElementById("analysisResult").innerHTML = `
            <div class="prediction ${cls}">
                <h2>${d.prediction === "R" ? "🔴 RESISTANT (R)" : "🟢 SUSCEPTIBLE (S)"}</h2>
                <p>${d.reason}</p>
            </div>
            <div class="card" style="margin-top:15px;">
                <p><b>Detected Determinants:</b> ${genes}</p>
                <p><b>Comparison Status:</b> ${d.status}</p>
            </div>
        `;
    } catch(e) {
        document.getElementById("analysisResult").innerHTML = `<div style="color:red;">Failed: ${e}</div>`;
    } finally {
        document.getElementById("loader").style.display = "none";
    }
}

loadSummary();
</script>
</body>
</html>"""

@app.get("/", response_class=HTMLResponse)
def index(): return HTML_PAGE

@app.get("/api/summary")
def api_summary():
    out = dict(metrics)
    for k in ["accuracy", "sensitivity", "specificity", "precision"]:
        out[k] = round(out[k] * 100, 2)
    return out

@app.get("/api/dataset")
def api_dataset(search: str = ""):
    if validation_df.empty: return []
    df = validation_df
    if search:
        df = df[df["BV_BRC_genome_id"].astype(str).str.lower().str.contains(search.lower(), na=False)]
    cols = ["BV_BRC_genome_id", "Meropenem_AST", "GenCon_Genes_Detected", "GenCon_Prediction", "Status"]
    return [row_to_json(r) for _, r in df.head(361)[[c for c in cols if c in df.columns]].iterrows()]

@app.post("/api/analyze-fasta")
async def api_analyze_fasta(file: UploadFile = File(...), ast: str = Form("")):
    tdir = tempfile.mkdtemp(prefix="gencon_")
    fpath = os.path.join(tdir, file.filename or "genome.fasta")
    try:
        content = await file.read()
        with open(fpath, "wb") as f:
            f.write(content)
        return run_fasta_analysis(fpath, ast)
    finally:
        shutil.rmtree(tdir, ignore_errors=True)
