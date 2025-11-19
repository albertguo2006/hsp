import pandas as pd
from bioservices import UniProt
import requests
import logging
import io

# Import configuration
try:
    import config
except ImportError:
    print("Warning: 'config.py' not found. Using default settings.")
    class config:
        BIOGRID_ACCESS_KEY = "YOUR_ACCESS_KEY_HERE"

# ==========================================
# CONFIGURATION
# ==========================================
CHAPERONE_LIST = [
    "DJC24_HUMAN", "DJC27_HUMAN", "DNJ5B_HUMAN", "DNJ5G_HUMAN",
    "DNJA1_HUMAN", "DNJA2_HUMAN", "DNJA3_HUMAN", "DNJA4_HUMAN",
    "DNJB1_HUMAN", "DNJB2_HUMAN", "DNJB3_HUMAN", "DNJB4_HUMAN",
    "DNJB5_HUMAN", "DNJB6_HUMAN", "DNJB7_HUMAN", "DNJB8_HUMAN",
    "DNJC2_HUMAN", "DNJC5_HUMAN", "DNJC7_HUMAN", "DNJC8_HUMAN",
    "DNJC9_HUMAN", "SACS_HUMAN"
]

USE_BIOGRID = False 
BIOGRID_ACCESS_KEY = config.BIOGRID_ACCESS_KEY

INTERACTIONS_FILE = "chaperone_interactions.csv"
DOMAINS_FILE = "target_domains.csv"
MASTER_FILE = "chaperone_domain_analysis_master.csv"

# Silence Bioservices logging
logging.getLogger("bioservices").setLevel(logging.ERROR)
u = UniProt()

def get_uniprot_accessions(gene_names):
    print(f"Mapping {len(gene_names)} chaperone names to UniProt IDs...")
    results = u.mapping("UniProtKB_AC-ID", "UniProtKB", ",".join(gene_names))
    mapping_dict = {}
    
    if "results" in results:
        for result in results["results"]:
            input_name = result["from"]
            val = result["to"]
            
            # FIX: Handle case where 'to' is a dictionary object
            if isinstance(val, dict):
                # Try to find primary accession
                accession = val.get("primaryAccession", str(val))
            else:
                accession = str(val)
                
            mapping_dict[input_name] = accession
            
    print(f"Successfully mapped {len(mapping_dict)} chaperones.")
    return mapping_dict

def get_intact_interactions(chaperone_acc):
    """
    Queries IntAct via PSICQUIC web service for protein interactions.
    Documentation: https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/
    """
    targets = set()
    
    # Correct PSICQUIC endpoint for IntAct
    intact_url = f"https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{chaperone_acc}?format=tab25"
    
    try:
        r = requests.get(intact_url, timeout=30)
        
        if r.status_code != 200:
            print(f"  Error: IntAct API returned status {r.status_code}")
            return []
            
        lines = r.text.splitlines()
        
        for line in lines:
            if not line or line.startswith('#'): continue
                
            cols = line.split('\t')
            if len(cols) < 2: continue
            
            # Parse Columns 0 and 1 (ID A and ID B)
            raw_id_a = cols[0]
            raw_id_b = cols[1]
            
            def parse_id(raw):
                if "uniprotkb:" in raw:
                    return raw.split(":")[1].split("-")[0]
                return None

            id_a = parse_id(raw_id_a)
            id_b = parse_id(raw_id_b)
            
            if id_a and id_a != chaperone_acc: targets.add(id_a)
            if id_b and id_b != chaperone_acc: targets.add(id_b)
            
    except Exception as e:
        print(f"  Error querying IntAct: {e}")
        
    return list(targets)

def get_biogrid_interactions(chaperone_acc):
    if not USE_BIOGRID or not BIOGRID_ACCESS_KEY or BIOGRID_ACCESS_KEY == "YOUR_ACCESS_KEY_HERE":
        return []
    targets = set()
    url = "https://webservice.thebiogrid.org/interactions/"
    params = {
        "accessKey": BIOGRID_ACCESS_KEY,
        "format": "json",
        "searchNames": "true",
        "geneList": chaperone_acc,
        "includeInteractors": "true",
        "taxId": 9606 
    }
    try:
        r = requests.get(url, params=params)
        if r.status_code == 200:
            data = r.json()
            # BioGRID logic placeholder
            pass 
    except Exception as e:
        print(f"  BioGRID Error: {e}")
    return list(targets)

def get_domains_batch(accession_list):
    print(f"Retrieving domains for {len(accession_list)} unique targets...")
    domain_data = []
    
    chunk_size = 20 
    for i in range(0, len(accession_list), chunk_size):
        chunk = accession_list[i:i+chunk_size]
        chunk = [x for x in chunk if len(x) < 15] 
        if not chunk: continue

        query = " OR ".join([f"accession:{acc}" for acc in chunk])
        
        try:
            df_res = u.get_df(query, columns="id,protein_name,xref_pfam")
            
            if df_res is not None and not df_res.empty:
                for _, row in df_res.iterrows():
                    target_id = row.get("Entry")
                    # SAFETY FIX: Force string conversion
                    target_name = str(row.get("Protein names", "Unknown"))
                    pfam_entry = str(row.get("Pfam", ""))
                    
                    if not pfam_entry or pfam_entry.lower() == "nan":
                        continue
                    
                    pfams = [x for x in pfam_entry.split(";") if x.strip()]
                    
                    for pfam_id in pfams:
                        domain_data.append({
                            "Target_ID": target_id,
                            "Target_Protein_Name": target_name.split("(")[0].strip(),
                            "Domain_ID": pfam_id
                        })
        except Exception as e:
            print(f"  Error fetching batch {i}-{i+chunk_size}: {e}")
            
    return pd.DataFrame(domain_data)

# ==========================================
# MAIN EXECUTION
# ==========================================
def main():
    print("--- Starting Chaperone Domain Miner ---")
    
    chap_map = get_uniprot_accessions(CHAPERONE_LIST)
    interaction_data = []
    all_targets = set()
    
    print("\nMining interactions...")
    for name, acc in chap_map.items():
        print(f"  Processing {name} ({acc})...")
        
        # Ensure 'acc' is a string before passing
        if not isinstance(acc, str):
            print(f"    Skipping {name}: Invalid accession format {acc}")
            continue
            
        partners = get_intact_interactions(acc)
        print(f"    Found {len(partners)} partners in IntAct")
        
        if USE_BIOGRID:
             bg_partners = get_biogrid_interactions(acc)
             if bg_partners:
                partners.extend(bg_partners)
            
        for partner in partners:
            if partner != acc: 
                interaction_data.append({
                    "Chaperone_Name": name,
                    "Chaperone_ID": acc,
                    "Target_ID": partner,
                    "Source": "IntAct"
                })
                all_targets.add(partner)

    df_interactions = pd.DataFrame(interaction_data)
    if df_interactions.empty:
        print("\nNo interactions found. Exiting.")
        return

    df_interactions.to_csv(INTERACTIONS_FILE, index=False)
    print(f"\nSaved {len(df_interactions)} interactions to {INTERACTIONS_FILE}")

    unique_targets = list(all_targets)
    df_domains = get_domains_batch(unique_targets)
    
    if df_domains.empty:
        print("No domain information found for targets.")
        return

    df_domains.to_csv(DOMAINS_FILE, index=False)
    print(f"Saved domain info to {DOMAINS_FILE}")
    
    print("\nMerging data...")
    df_master = pd.merge(df_interactions, df_domains, on="Target_ID", how="inner")
    df_master["Interaction_Label"] = 1
    
    cols = ["Chaperone_Name", "Chaperone_ID", "Target_ID", "Target_Protein_Name", "Domain_ID", "Interaction_Label", "Source"]
    cols = [c for c in cols if c in df_master.columns]
    
    df_master = df_master[cols]
    df_master.to_csv(MASTER_FILE, index=False)
    
    print(f"--- SUCCESS ---")
    print(f"Master dataset created: {MASTER_FILE}")
    print(f"Total Rows: {len(df_master)}")

if __name__ == "__main__":
    main()