"""Turn the openFDA counts into PRR / ROR disproportionality signals, charts and a dashboard."""
import json, math, sqlite3
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.express as px

meta = json.load(open("data/meta.json")); N = meta["all_reports"]
con = sqlite3.connect("data/faers.db")
pd.read_csv("data/drug_totals.csv").to_sql("drug_totals", con, if_exists="replace", index=False)
pd.read_csv("data/drug_reaction_counts.csv").to_sql("drug_reaction", con, if_exists="replace", index=False)
pd.read_csv("data/reaction_totals.csv").to_sql("reaction_totals", con, if_exists="replace", index=False)
con.execute("DROP TABLE IF EXISTS meta"); con.execute("CREATE TABLE meta(n_all INTEGER)"); con.execute("INSERT INTO meta VALUES (?)", (N,))

# 2x2 table pieces in SQL: a = drug+reaction, drug total, reaction total, all reports
SQL = """
SELECT dr.drug, dr.reaction, dr.a,
       dt.drug_reports - dr.a            AS b,
       rt.reaction_reports - dr.a        AS c,
       m.n_all - dt.drug_reports - rt.reaction_reports + dr.a AS d,
       dt.drug_reports, rt.reaction_reports
FROM drug_reaction dr
JOIN drug_totals dt   ON dt.drug = dr.drug
JOIN reaction_totals rt ON rt.reaction = dr.reaction
CROSS JOIN meta m
"""
df = pd.read_sql(SQL, con)
df["PRR"] = (df.a / (df.a + df.b)) / (df.c / (df.c + df.d))
df["ROR"] = (df.a * df.d) / (df.b * df.c)
se_prr = np.sqrt(1/df.a - 1/(df.a+df.b) + 1/df.c - 1/(df.c+df.d))
se_ror = np.sqrt(1/df.a + 1/df.b + 1/df.c + 1/df.d)
df["PRR_low"] = np.exp(np.log(df.PRR) - 1.96*se_prr)
df["ROR_low"] = np.exp(np.log(df.ROR) - 1.96*se_ror)
# Evans criteria: PRR >= 2, at least 3 reports, chi-square >= 4 (here: lower 95% bound of PRR > 1 as a stand-in)
df["signal"] = (df.PRR >= 2) & (df.a >= 3) & (df.PRR_low > 1)
df["share_of_drug_reports"] = df.a / df.drug_reports
df = df.sort_values(["drug", "PRR"], ascending=[True, False])
df.round(3).to_csv("data/signals_all.csv", index=False)
sig = df[df.signal]
sig.round(3).to_csv("data/signals_flagged.csv", index=False)
print(df.groupby("drug").signal.sum())

# Figure 1: top 8 signals per drug (ranked by PRR, min 200 reports)
top = (sig[sig.a >= 200].sort_values("PRR", ascending=False).groupby("drug").head(5))
fig, ax = plt.subplots(figsize=(10, 12))
labels = top.drug.str.title() + " | " + top.reaction.str.lower()
order = top.sort_values(["drug", "PRR"]).reset_index(drop=True)
colors = dict(zip(sorted(top.drug.unique()), plt.cm.tab20.colors))
ax.barh(range(len(order)), order.PRR, color=[colors[d] for d in order.drug])
ax.set_yticks(range(len(order)))
ax.set_yticklabels(order.drug.str.title() + ": " + order.reaction.str.lower(), fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("PRR (log scale)")
ax.set_title("Strongest reaction signals per drug (PRR, at least 200 reports)\nopenFDA FAERS, reports received 2023-2025")
plt.tight_layout(); plt.savefig("figures/top_signals_per_drug.png", dpi=150); plt.close()

# Figure 2: heatmap of log PRR for the most-reported reactions across drugs
common = df.groupby("reaction").a.sum().sort_values(ascending=False).head(25).index
piv = df[df.reaction.isin(common)].pivot(index="reaction", columns="drug", values="PRR")
fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(np.log2(piv.values), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns.str.title(), rotation=45, ha="right")
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index.str.lower(), fontsize=8)
plt.colorbar(im, label="log2(PRR): red = reported more than usual, blue = less")
ax.set_title("Common reactions across 12 drugs (blank = not in that drug's top 60)")
plt.tight_layout(); plt.savefig("figures/reaction_heatmap.png", dpi=150); plt.close()

# Figure 3: share of reports that mention drug-ineffective style terms vs GI terms, per drug
gi = ["NAUSEA", "VOMITING", "DIARRHOEA", "DECREASED APPETITE", "CONSTIPATION", "ABDOMINAL PAIN UPPER"]
g = df[df.reaction.isin(gi)].groupby("drug").share_of_drug_reports.sum().sort_values()
fig, ax = plt.subplots(figsize=(8, 5))
ax.barh(g.index.str.title(), g.values * 100, color="#2a7f9e")
ax.set_xlabel("Gastrointestinal terms per 100 reports (sum of 6 terms, a report can count more than once)")
ax.set_title("GI burden in the reports, by drug")
plt.tight_layout(); plt.savefig("figures/gi_share_by_drug.png", dpi=150); plt.close()

# Interactive dashboard
d2 = df[df.a >= 100].copy(); d2["PRR_plot"] = d2.PRR.clip(upper=100)
f = px.scatter(d2, x="a", y="PRR_plot", color="drug", hover_name="reaction", log_x=True, log_y=True,
               hover_data={"PRR": ":.2f", "ROR": ":.2f", "a": True, "PRR_plot": False, "signal": True},
               labels={"a": "reports with this drug and reaction", "PRR_plot": "PRR (capped at 100)"},
               title="Drug safety signals in FAERS 2023-2025 (hover for the reaction)")
f.add_hline(y=2, line_dash="dash", annotation_text="PRR = 2")
f.write_html("dashboard.html", include_plotlyjs="cdn")
print("flagged", len(sig), "of", len(df))
