import pandas as pd
import os
import sys

# Configuration
MASTER_FILE = "chaperone_domain_analysis_master.csv"

def analyze_preferences():
    # 1. Check if data exists
    if not os.path.exists(MASTER_FILE):
        print(f"Error: '{MASTER_FILE}' not found.")
        print("Please run 'chaperone_miner.py' first to generate the data.")
        return

    # 2. Load the data
    print(f"Loading {MASTER_FILE}...")
    try:
        df = pd.read_csv(MASTER_FILE)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # 3. Basic Stats
    total_interactions = len(df)
    unique_chaperones = df['Chaperone_Name'].unique()
    print(f"Loaded {total_interactions} domain interactions for {len(unique_chaperones)} chaperones.\n")
    
    # 4. Analyze each chaperone
    for chap in unique_chaperones:
        print(f"=========================================")
        print(f"  PREFERENCES FOR: {chap}")
        print(f"=========================================")
        
        # Filter data for just this chaperone
        subset = df[df['Chaperone_Name'] == chap]
        total_targets = len(subset)
        
        if total_targets == 0:
            print("  No targets found.\n")
            continue
            
        # Count domain frequencies
        # We group by Domain_ID and count how many times it appears
        domain_counts = subset['Domain_ID'].value_counts()
        
        print(f"  Total Domains Found: {total_targets}")
        print(f"  Top 15 Most Frequent Domains:")
        print(f"  ---------------------------------------------")
        print(f"  {'Domain ID':<12} | {'Count':<6} | {'% of Clients':<12}")
        print(f"  ---------------------------------------------")
        
        for domain, count in domain_counts.head(15).items():
            percent = (count / total_targets) * 100
            # Basic interpretation helper (You can expand this dictionary manually if specific IDs are common)
            # e.g. PF00069 is usually Protein Kinase
            print(f"  {domain:<12} | {count:<6} | {percent:.1f}%")
        
        print("\n")

if __name__ == "__main__":
    analyze_preferences()