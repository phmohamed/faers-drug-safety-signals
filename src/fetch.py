"""Pull adverse event counts from the openFDA FAERS endpoint and save them as CSV.
Window: reports received 2023-01-01 to 2025-12-31. Needs internet, no API key."""
import json, time, urllib.parse, urllib.request, csv, os
BASE = "https://api.fda.gov/drug/event.json"
WINDOW = "receivedate:[20230101+TO+20251231]"
DRUGS = ["semaglutide","tirzepatide","metformin","empagliflozin","dapagliflozin","sitagliptin",
         "atorvastatin","apixaban","adalimumab","pembrolizumab","sertraline","ibuprofen"]
TOP_N = 60
def get(search, count=None, limit=1):
    url = f"{BASE}?search={search}"
    if count: url += f"&count={count}&limit={limit}"
    else: url += "&limit=1"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url.replace(" ","%20"), timeout=40) as r:
                return json.load(r)
        except Exception as e:
            time.sleep(2 + attempt * 2)
    return None
def total(search):
    d = get(search)
    return d["meta"]["results"]["total"] if d else None
os.makedirs("data", exist_ok=True)
N = total(WINDOW)
print("all reports in window:", N)
drug_rows, reaction_rows = [], []
for dname in DRUGS:
    s = f'{WINDOW}+AND+patient.drug.openfda.generic_name:"{dname}"'
    nd = total(s)
    d = get(s, "patient.reaction.reactionmeddrapt.exact", TOP_N)
    print(dname, nd)
    for r in d["results"]:
        reaction_rows.append((dname, r["term"], r["count"]))
    drug_rows.append((dname, nd))
    time.sleep(0.3)
terms = sorted({t for _, t, _ in reaction_rows})
CACHE = "data/.reaction_cache.json"
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
t0 = time.time()
for t in terms:
    if t in cache: continue
    if time.time() - t0 > 85:
        json.dump(cache, open(CACHE, "w")); print("partial, rerun"); raise SystemExit
    s = f'{WINDOW}+AND+patient.reaction.reactionmeddrapt.exact:"{urllib.parse.quote(t)}"'
    cache[t] = total(s)
json.dump(cache, open(CACHE, "w"))
tot_rows = [(t, cache[t]) for t in terms]
with open("data/drug_totals.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["drug","drug_reports"]); w.writerows(drug_rows)
with open("data/drug_reaction_counts.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["drug","reaction","a"]); w.writerows(reaction_rows)
with open("data/reaction_totals.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["reaction","reaction_reports"]); w.writerows(tot_rows)
json.dump({"all_reports": N, "window": "2023-01-01 to 2025-12-31", "source": "openFDA drug/event"}, open("data/meta.json","w"))
print("done", len(terms), "reaction terms")
