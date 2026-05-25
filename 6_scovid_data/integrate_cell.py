import pandas as pd
import os
import glob

# Comprehensive list of cell types found in your directory tree
TARGET_CELL_TYPES = [
    'cell_B', 'cell_Endothelial', 'cell_Epithelial', 'cell_Fibroblasts',
    'cell_Macrophages', 'cell_Monocytes', 'cell_Neuronal', 'cell_SmoothMuscle',
    'cell_T', 'cell_Unknown', 'cell_DC', 'cell_Mast', 'cell_Neutrophils', 
    'cell_RBC', 'cell_Glial', 'cell_Brain', 'cell_Monoctyes' 
]

def integrate_data():
    """
    Integrates CSVs from lung, airway, and brain info folders.
    Tracks source tissue, GEO dataset, and standardizes cell types.
    """
    all_data = []
    info_folders = ['lung_info', 'airway_info', 'brain_info']
    
    print("=== STARTING CHROMIUM DATA INTEGRATION ===\n")

    for info_dir in info_folders:
        if not os.path.exists(info_dir):
            print(f"⚠️  Skipping: {info_dir} (Directory not found)")
            continue
            
        folder_csv_count = 0
        print(f"📁 Processing Tissue: {info_dir}")
        
        # Identify GSE subdirectories
        gse_subdirs = [d for d in os.listdir(info_dir) if os.path.isdir(os.path.join(info_dir, d)) and d.startswith('GSE')]
        
        for gse_dir in gse_subdirs:
            current_path = os.path.join(info_dir, gse_dir)
            csv_files = glob.glob(os.path.join(current_path, "*.csv"))
            
            print(f"   └── 📂 {gse_dir}: Found {len(csv_files)} CSVs")
            
            for csv_file in csv_files:
                filename = os.path.basename(csv_file)
                
                try:
                    df = pd.read_csv(csv_file)
                    
                    # Logic to extract cell type from filename
                    # Example: GSE171524_cell_DC.csv -> cell_DC
                    filename_no_ext = filename.replace('.csv', '')
                    file_key = filename_no_ext.replace(f"{gse_dir}_", "")
                    
                    if 'total' in filename.lower():
                        if 'cell_type' not in df.columns:
                            df['cell_type'] = 'total'
                    else:
                        # Ensure 'cell_' prefix exists
                        cell_type = file_key if file_key.startswith('cell_') else f"cell_{file_key}"
                        df['cell_type'] = cell_type
                    
                    # Add metadata for integration
                    df['source_tissue'] = info_dir.replace('_info', '')
                    df['GEO_ID'] = gse_dir
                    df['origin_filename'] = filename
                    
                    # Clean column names (strip whitespace and quotes)
                    df.columns = [col.replace('"', '').strip() for col in df.columns]
                    
                    all_data.append(df)
                    folder_csv_count += 1
                    
                except Exception as e:
                    print(f"   ❌ Error processing {filename}: {str(e)}")
        
        print(f"✅ Finished {info_dir}: {folder_csv_count} files processed.\n")

    if all_data:
        # Merge everything
        combined_df = pd.concat(all_data, ignore_index=True)
        
        # Save Master File
        master_filename = 'smerged_covid.csv'
        combined_df.to_csv(master_filename, index=False)
        
        print("-" * 40)
        print(f"🚀 SUCCESS: Integrated {len(combined_df)} rows.")
        print(f"💾 Master file: {master_filename}")
        
        # Export grouped files by cell type
        print("\nCreating grouped files by cell type...")
        # Use unique values found to ensure even types not in our list get saved
        found_types = combined_df['cell_type'].unique()
        
        for c_type in found_types:
            subset = combined_df[combined_df['cell_type'] == c_type]
            output_path = f'grouped_{c_type}.csv'
            subset.to_csv(output_path, index=False)
            print(f"   📦 {output_path} ({len(subset)} rows)")

        print("\n=== FINAL SUMMARY BY TISSUE ===")
        print(combined_df['source_tissue'].value_counts())
        
        return combined_df
    else:
        print("🛑 Error: No data found. Check your root directory path.")
        return None

if __name__ == "__main__":
    integrate_data()