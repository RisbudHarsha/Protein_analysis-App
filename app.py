import streamlit as st
import pandas as pd
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -----------------------------
# PAGE CONFIG
# -----------------------------
st.set_page_config(page_title="Protein Analysis Suite", layout="wide")

# -----------------------------
# UI STYLING
# -----------------------------
st.markdown("""
<style>
.stApp {
    background: linear-gradient(rgba(13,17,23,0.85), rgba(13,17,23,0.95)),
                url("https://images.unsplash.com/photo-1532187863486-abf9dbad1b69");
    background-size: cover;
    background-position: center;
    color: #e6edf3;
}

.card {
    background-color: rgba(28, 33, 40, 0.85);
    padding: 20px;
    border-radius: 14px;
    margin-bottom: 15px;
    border: 1px solid #30363d;
    backdrop-filter: blur(6px);
}

h1 {color: #58a6ff;}
h2 {color: #79c0ff;}
h3 {color: #a5d6ff;}

.stButton>button {
    background-color: #238636;
    color: white;
    border-radius: 10px;
    font-weight: 600;
}
</style>
""", unsafe_allow_html=True)

# -----------------------------
# FUNCTIONS
# -----------------------------
def clean_sequence(seq):
    return seq.replace(" ", "").replace("\n", "").upper().replace("U", "C").replace("O", "K")

def fetch_sequence(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    r = requests.get(url)
    if r.status_code == 200:
        return "".join(r.text.split("\n")[1:])
    return None

def fetch_isoforms(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    r = requests.get(url)

    isoforms = []
    if r.status_code == 200:
        data = r.json()
        for c in data.get("comments", []):
            if c["commentType"] == "ALTERNATIVE PRODUCTS":
                for iso in c.get("isoforms", []):
                    iso_id = iso.get("isoformIds", [None])[0]
                    if iso_id:
                        fasta_url = f"https://rest.uniprot.org/uniprotkb/{iso_id}.fasta"
                        res = requests.get(fasta_url)
                        if res.status_code == 200:
                            seq = "".join(res.text.split("\n")[1:])
                            isoforms.append((iso_id, seq))
    return isoforms

def analyze(seq):
    p = ProteinAnalysis(clean_sequence(seq))
    return {
        "Molecular Weight": round(p.molecular_weight(), 2),
        "Isoelectric Point": round(p.isoelectric_point(), 2),
        "Instability Index": round(p.instability_index(), 2),
        "GRAVY": round(p.gravy(), 2)
    }

def functional_interpretation(res):
    insights = []
    insights.append("This protein seems stable under typical conditions." if res["Instability Index"] < 40 
                    else "This protein may be less stable and could degrade more easily.")
    insights.append("It leans towards hydrophobic behavior." if res["GRAVY"] > 0 
                    else "It likely interacts well with water (hydrophilic).")
    insights.append("Overall, it behaves like a basic protein." if res["Isoelectric Point"] > 7 
                    else "Overall, it behaves like an acidic protein.")
    return insights

def secondary_structure(seq):
    p = ProteinAnalysis(clean_sequence(seq))
    h, t, s = p.secondary_structure_fraction()
    return {"Helix": h, "Turn": t, "Sheet": s}

def compare_isoforms(isoforms):
    data = []
    for name, seq in isoforms:
        res = analyze(seq)
        res["Isoform"] = name
        data.append(res)
    return pd.DataFrame(data).set_index("Isoform")

def generate_report(res, insights):
    text = "Protein Analysis Summary\n\n"
    for k, v in res.items():
        text += f"{k}: {v}\n"
    text += "\nInterpretation:\n"
    for i in insights:
        text += f"- {i}\n"
    return text

# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.title("🧬 Protein Analysis Suite")

section = st.sidebar.radio("Where would you like to go?", ["Tool", "About", "Team"])

st.sidebar.markdown("---")
st.sidebar.subheader("Quick Start")
st.sidebar.write("""
1. Pick an analysis type  
2. Provide your sequence or UniProt ID  
3. Run the analysis  
4. Explore the results  
""")

# =============================
# TOOL
# =============================
if section == "Tool":

    tool = st.sidebar.selectbox(
        "What would you like to explore?",
        ["Protein Analysis", "Isoform Analysis"]
    )

    st.title("🧬 Protein Analysis Dashboard")
    st.markdown("Take a closer look at your protein and uncover its key biochemical features.")

    # -------------------------
    # PROTEIN ANALYSIS
    # -------------------------
    if tool == "Protein Analysis":

        mode = st.radio("How would you like to begin?", ["Paste Sequence", "Use UniProt ID"])

        if mode == "Paste Sequence":
            seq = st.text_area("You can paste your protein sequence here")
        else:
            uid = st.text_input("Enter a UniProt ID")
            seq = fetch_sequence(uid) if uid else None

        if st.button("Run Analysis"):

            if not seq:
                st.error("I couldn’t find a valid sequence. Please double-check your input and try again.")
            else:
                with st.spinner("Working on your sequence... just a moment"):
                    res = analyze(seq)
                    insights = functional_interpretation(res)

                st.success("Analysis complete!")

                # RESULTS
                st.markdown('<div class="card">', unsafe_allow_html=True)
                st.subheader("Key Properties")
                st.write("Here’s a quick overview of the main physicochemical properties:")
                st.dataframe(pd.DataFrame(res.items(), columns=["Property", "Value"]))
                st.markdown('</div>', unsafe_allow_html=True)

                # SECONDARY STRUCTURE
                st.markdown('<div class="card">', unsafe_allow_html=True)
                st.subheader("Secondary Structure Snapshot")
                st.write("This gives an estimate of how the protein folds (helix, turn, sheet):")
                st.bar_chart(pd.Series(secondary_structure(seq)))
                st.markdown('</div>', unsafe_allow_html=True)

                # INTERPRETATION
                st.markdown('<div class="card">', unsafe_allow_html=True)
                st.subheader("What can we take from this?")
                for i in insights:
                    st.write(f"• {i}")
                st.markdown('</div>', unsafe_allow_html=True)

                st.markdown("💡 Tip: You can cross-check these results on ExPASy ProtParam.")

                report = generate_report(res, insights)
                st.download_button("Download Summary", report)

    # -------------------------
    # ISOFORM ANALYSIS
    # -------------------------
    else:

        uid = st.text_input("Enter a UniProt ID to compare isoforms")

        if st.button("Compare Isoforms"):

            with st.spinner("Fetching isoform data..."):
                isoforms = fetch_isoforms(uid)

            if len(isoforms) < 2:
                st.warning("I couldn’t find enough isoforms for comparison.")
            else:
                df = compare_isoforms(isoforms)

                st.success("Isoforms loaded successfully!")

                st.subheader("Comparison Table")
                st.dataframe(df)

                st.subheader("Visual Comparison")

                col1, col2 = st.columns(2)

                with col1:
                    st.write("Molecular Weight")
                    st.bar_chart(df["Molecular Weight"])

                    st.write("Isoelectric Point")
                    st.bar_chart(df["Isoelectric Point"])

                with col2:
                    st.write("Instability Index")
                    st.bar_chart(df["Instability Index"])

                    st.write("GRAVY")
                    st.bar_chart(df["GRAVY"])

# =============================
# ABOUT
# =============================
elif section == "About":
    st.title("About this tool")

    st.write("""
This application is designed to help you quickly analyze protein sequences and understand their basic properties.

### 🔬 What this tool does
• Calculates key physicochemical properties (molecular weight, pI, stability, hydrophobicity)  
• Predicts secondary structure composition (helix, turn, sheet)  
• Provides simple interpretation of results  
• Allows comparison of protein isoforms  

### ⚙️ How it works
• Uses BioPython (ProtParam) for calculations  
• Fetches sequence data from UniProt  
• Displays results in tables and charts  

### 🎯 Who can use it
• Students learning bioinformatics  
• Beginners exploring protein analysis  
• For quick academic or research insights  

### 💡 Why this tool
• Simple and easy to use   
• Combines analysis + visualization in one place  
""")
# =============================
# TEAM
# =============================
elif section == "Team":
    st.title("Team")
    st.write("""
Built with curiosity and learning by:

**Akash Roy**  
**Snehal Yadav**  
**Harshada Risbud**  

""")