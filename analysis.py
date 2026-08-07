import pandas as pd
from scipy import stats
import sys
import os
import re
import time
import datetime
from tqdm import tqdm
import warnings
import logging
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import threading
from statsmodels.stats.multitest import multipletests
import itertools
import pyperclip
import numpy as np

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("\n[System Tip] adjustText is not installed (pip install adjustText), volcano plot labels may overlap.")

# ================= [v81.1] Global Environment and Plot Settings =================
IMG_EXT = "png"  
IMG_DPI = 300    
FDR_ENABLED = True
DIST_PLOT_STYLE = "3"
SESSION_TIME = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
SESSION_FOLDER_NAME = f"GeneAnalysis_{SESSION_TIME}"

def initialize_environment():
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message=".*Glyph.*")
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=stats.ConstantInputWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Precision loss.*")
    
    # Suppress ConvergenceWarning for complete separation in survival analysis
    try:
        from lifelines.utils import ConvergenceWarning
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
    except:
        pass

    plt.rcdefaults()
    sns.set(style="ticks")
    fonts = ['Microsoft JhengHei', 'SimHei', 'Arial', 'Helvetica', 'sans-serif']
    plt.rcParams.update({
        'font.sans-serif': fonts,
        'axes.unicode_minus': False,
        'figure.dpi': IMG_DPI,
        'axes.titlesize': 22,
        'axes.labelsize': 20,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 16,
        'legend.title_fontsize': 18,
        'font.size': 16
    })

initialize_environment()
# ==========================================================

# ================= 1. Spinner Animation and Helper Inputs =================
class Spinner:
    def __init__(self, message="Loading...", delay=0.1):
        self.spinner = itertools.cycle(['|', '/', '-', '\\'])
        self.delay = delay
        self.busy = False
        self.spinner_visible = False
        self.message = message

    def run(self):
        while self.busy:
            sys.stdout.write(f'\r{self.message} ' + next(self.spinner))
            sys.stdout.flush()
            time.sleep(self.delay)

    def __enter__(self):
        self.busy = True
        self.thread = threading.Thread(target=self.run)
        self.thread.start()

    def __exit__(self, exception, value, tb):
        self.busy = False
        self.thread.join()
        if exception:
            sys.stdout.write(f'\r{self.message} [Failed] ❌\n')
        else:
            sys.stdout.write(f'\r{self.message} [Done] ✅  \n')
        sys.stdout.flush()

def get_input_list(prompt_text, allow_none=False):
    print(prompt_text)
    if allow_none:
        print("   [Tip] Press Enter to read from clipboard (or type N to skip); you can also type directly (comma-separated)")
        print("   [Tip] You can copy a list from Excel/web/text file, then press Enter to read from clipboard (or type N to skip)")
    else:
        print("   [Tip] Press Enter to read from clipboard; you can also type directly (comma-separated)")
        print("   [Tip] You can copy a list from Excel/web/text file, then press Enter to read from clipboard; you can also type directly (comma-separated)")
    
    user_input = input(" ➤ ").strip()
    
    if allow_none and user_input.upper() in ['N', 'NONE', 'NO']:
        return []
        
    if user_input == '' or user_input.upper() == 'C':
        try:
            content = pyperclip.paste()
            if content.strip():
                print(f" 📋 [System] Content read from clipboard: {content[:30]}...")
                user_input = content
            else:
                if allow_none: return []
                print(" ⚠️ [System] Clipboard is empty, please re-enter.")
                return []
        except Exception:
            if allow_none: return []
            print(" ⚠️ [System] Failed to read from clipboard, please paste manually.")
            return []
            
    clean_text = user_input.replace('\u3000', ' ')
    tokens = re.split(r'[,;\t\n|]+', clean_text)
    
    final_list = [t.strip() for t in tokens if t.strip()]
    if final_list:
        print(f" -> Successfully identified {len(final_list)} items")
    return final_list

def get_default_input(prompt_text, default_val):
    user_in = input(f"{prompt_text} [Default: '{default_val}']: ").strip()
    return user_in if user_in else default_val

def play_beep(times=3, delay=0.3):
    """Play system beep sound multiple times"""
    for _ in range(times):
        sys.stdout.write('\a')
        sys.stdout.flush()
        time.sleep(delay)

def normalize_name(name):
    return str(name).strip().upper().replace('_', ' ')

def ask_transform_method():
    print("\nPlease select data transformation/calculation method (will be applied to all genes and targets in this session):")
    print("  (1) Raw data (default, preserves true linear distance)")
    print("  (2) Winsorization (clips outliers to -3 ~ +3, prevents chart compression)")
    print("  (3) Percentile Rank (absolutely uniform, necessary for non-parametric stats)")
    ans = input(">> ").strip()
    if ans not in ['1', '2', '3']: ans = '1'
    
    if ans == '2': print("\n🛡️ [System] Winsorization (-3 ~ +3) has been applied globally")
    elif ans == '3': print("\n📊 [System] Percentile Rank transformation (Non-parametric) has been applied globally")
    
    return ans

def transform_series(s, method):
    if s is None: return None
    if method == '2': 
        s_num = pd.to_numeric(s, errors='coerce')
        s_z = (s_num - s_num.mean()) / (s_num.std() + 1e-9)
        return s_z.clip(lower=-3, upper=3)
    if method == '3': return s.rank(pct=True)
    return s

def transform_df(df, method):
    if df is None: return None
    if method == '2': 
        df_num = df.apply(pd.to_numeric, errors='coerce')
        df_z = (df_num - df_num.mean()) / (df_num.std() + 1e-9)
        return df_z.clip(lower=-3, upper=3)
    if method == '3': return df.rank(pct=True)
    return df

def get_p_stars(p):
    if pd.isna(p) or p is None:
        return ""
    if p < 0.0001: return "****"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."

def format_pval(p, use_e=True, with_stars=False):
    if pd.isna(p) or p is None:
        return "= N/A"
    
    # Determine the p-value string representation
    if p < 0.0001:
        p_str = "< 0.0001"
    elif use_e:
        p_str = f"= {p:.3e}"
    else:
        p_str = f"= {p:.3f}"
        
    if with_stars:
        stars = get_p_stars(p)
        if stars == "n.s.":
            return f"{p_str} (n.s.)"
        else:
            return f"{p_str} {stars}"
    else:
        return p_str

def ask_local_subtyping(df_clinical, default_subtype_col):
    local_subtype_col = default_subtype_col
    subtypes_to_run = None
    if df_clinical is not None:
        ans = input(f"\n➤ Current global grouping is '{default_subtype_col}'. Customize analysis groups for this run? (y/n) [Default: n]: ").strip().lower()
        if ans == 'y':
            selected_cols = []
            while True:
                level = len(selected_cols) + 1
                prompt = f"【Local Grouping {level}: Primary】Enter grouping column name (Press Enter to skip/finish): " if level == 1 else f"【Local Grouping {level}: Secondary】Sub-divide within {'_'.join(selected_cols)}? Enter column name (Press Enter to skip/finish): "
                c = input(prompt).strip()
                if not c:
                    break
                if c in df_clinical.columns:
                    selected_cols.append(c)
                else:
                    print(f"⚠️ Column '{c}' not found. Please re-enter.")

            if selected_cols:
                if len(selected_cols) == 1:
                    local_subtype_col = selected_cols[0]
                    print(f"-> 🎯 Local single-level grouping: {local_subtype_col}")
                else:
                    local_subtype_col = "_".join(selected_cols)
                    mask = df_clinical[selected_cols[0]].notna()
                    for col in selected_cols[1:]:
                        mask &= df_clinical[col].notna()
                    df_clinical.loc[mask, local_subtype_col] = df_clinical.loc[mask, selected_cols[0]].astype(str)
                    for col in selected_cols[1:]:
                        df_clinical.loc[mask, local_subtype_col] += "_" + df_clinical.loc[mask, col].astype(str)
                    print(f"-> 🎯 Local {len(selected_cols)}-level stratification enabled! Grouping by: {local_subtype_col}")

                avail_subs = sorted(list(df_clinical[local_subtype_col].dropna().unique()), key=str)
                print(f"Available subtypes: {', '.join(avail_subs)}")
                sub_sel = input("➤ Enter subtypes to run (comma-separated, 'all' for all, Enter for All_Patients only): ").strip()
                if sub_sel.lower() == 'all':
                    subtypes_to_run = None
                elif sub_sel:
                    subtypes_to_run = [s.strip() for s in re.split(r'[,;\t\n|]+', sub_sel) if s.strip() in avail_subs]
                    if not subtypes_to_run: subtypes_to_run = [None]
                else:
                    subtypes_to_run = [None]
            else:
                local_subtype_col = None
                subtypes_to_run = [None]
    return local_subtype_col, subtypes_to_run

def get_analysis_groups(df, subtype_col, subtypes_to_run=None):
    groups = []
    if subtypes_to_run is None:
        groups.append(('All_Patients', df))
        if subtype_col and subtype_col in df.columns:
            for st in sorted(list(df[subtype_col].dropna().unique()), key=str):
                groups.append((f"Sub_{st}", df[df[subtype_col] == st]))
    else:
        for st in subtypes_to_run:
            if st is None:
                groups.append(('All_Patients', df))
            elif subtype_col and subtype_col in df.columns:
                groups.append((f"Sub_{st}", df[df[subtype_col] == st]))
    return groups

def get_analysis_indices(df_clinical, target_index, subtype_col, subtypes_to_run=None):
    groups = []
    base_index = df_clinical.index if df_clinical is not None else target_index
    if subtypes_to_run is None:
        groups.append(('All_Patients', base_index))
        if df_clinical is not None and subtype_col and subtype_col in df_clinical.columns:
            for st in sorted(list(df_clinical[subtype_col].dropna().unique()), key=str):
                groups.append((f"Sub_{st}", df_clinical[df_clinical[subtype_col] == st].index))
    else:
        for st in subtypes_to_run:
            if st is None:
                groups.append(('All_Patients', base_index))
            elif df_clinical is not None and subtype_col and subtype_col in df_clinical.columns:
                groups.append((f"Sub_{st}", df_clinical[df_clinical[subtype_col] == st].index))
    return groups

def ask_for_palette(palette_type='diverging'):
    """
    Asks the user for a Matplotlib colormap/palette based on the chart type.
    palette_type: 'diverging', 'categorical', 'binary_hl'
    """
    if palette_type == 'diverging':
        print("\nSelect a diverging colormap for heatmaps/correlations:")
        print("  (1) vlag (Blue-White-Red, default)")
        print("  (2) coolwarm (Blue-White-Red, softer)")
        print("  (3) bwr (Blue-White-Red, stronger)")
        print("  (4) RdYlGn (Red-Yellow-Green)")
        print("  (5) PRGn (Purple-White-Green)")
        print("  (6) Custom (Enter a Matplotlib colormap name)")
        cmap_choice = input(">> ").strip() or '1'
        cmap_map = {'1': 'vlag', '2': 'coolwarm', '3': 'bwr', '4': 'RdYlGn', '5': 'PRGn'}
        if cmap_choice in cmap_map: return cmap_map[cmap_choice]
        if cmap_choice == '6': return input("Enter custom colormap name: ").strip() or 'vlag'
        return cmap_choice # Allow direct input of colormap name
    elif palette_type == 'categorical':
        print("\nSelect a categorical palette for group plots:")
        print("  (1) Set1 (Classic multi-color, default)")
        print("  (2) Pastel (Soft colors)")
        print("  (3) husl (Evenly spaced in HSL color space)")
        print("  (4) Custom (Enter a Matplotlib palette name)")
        pal_choice = input(">> ").strip() or '1'
        pal_map = {'1': 'Set1', '2': 'Pastel1', '3': 'husl'}
        if pal_choice in pal_map: return pal_map[pal_choice]
        if pal_choice == '4': return input("Enter custom palette name: ").strip() or 'Set1'
        return pal_choice # Allow direct input of palette name
    # Other types (e.g., binary_hl) can be added here
    return None

# ================= 2. Data Loading =================
def load_gene_file(file_path):
    try: return pd.read_csv(file_path, index_col=0, encoding='utf-8-sig').T
    except: return pd.read_csv(file_path, index_col=0, encoding='cp950').T

def load_pathway_file(file_path):
    try: df = pd.read_csv(file_path, index_col=0, encoding='utf-8-sig')
    except: df = pd.read_csv(file_path, index_col=0, encoding='cp950')
    return df.iloc[:, 1:].T 

def load_protein_file(file_path):
    try: df = pd.read_csv(file_path, index_col=0, encoding='utf-8-sig')
    except: df = pd.read_csv(file_path, index_col=0, encoding='cp950')
    if df.shape[1] > 0:
        try: pd.to_numeric(df.iloc[:, 0].dropna().iloc[0])
        except (ValueError, TypeError, IndexError):
            df = df.iloc[:, 1:]
    return df.T

def load_clinical_file(file_path):
    try: return pd.read_csv(file_path, index_col=0, encoding='utf-8-sig')
    except: return pd.read_csv(file_path, index_col=0, encoding='cp950')

def load_ext_data_basic(f1):
    if not os.path.exists(f1):
        print(f"\n[Error] File not found: {f1}")
        return None
    try: df = pd.read_csv(f1, encoding='utf-8-sig')
    except: df = pd.read_csv(f1, encoding='cp950')
    df.columns = df.columns.str.strip()
    possible_cols = ['Pathway', 'pathway', 'Term', 'term', 'TF', 'Gene', 'Symbol', 'ID', 'Name']
    gene_col = next((col for col in possible_cols if col in df.columns), df.columns[0])
    df.rename(columns={gene_col: 'Item_Name'}, inplace=True)
    
    if 'pvalue' in df.columns:
        df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
        min_positive_p = df.loc[df['pvalue'] > 0, 'pvalue'].min()
        if pd.isna(min_positive_p): min_positive_p = 1e-50
        
        # Add tiny random jitter to p=0 to prevent multiple extreme genes from overlapping on the Y-axis, which can cause adjustText to fail.
        zero_mask = df['pvalue'] == 0
        if zero_mask.sum() > 0:
            df.loc[zero_mask, 'pvalue'] = min_positive_p * np.random.uniform(0.05, 0.15, size=zero_mask.sum())
        
        df['mlog10p'] = -np.log10(df['pvalue'].astype(float))
    else:
        df['mlog10p'] = 0 
        
    return df

def load_ext_enrichment_data(f1, top_n, sort_col='mlog10p'):
    df = load_ext_data_basic(f1)
    if df is None: return None
    if sort_col not in df.columns and sort_col != 'mlog10p':
        print(f"[Warning] X-axis column '{sort_col}' not found in CSV, forcing to mlog10p")
        sort_col = 'mlog10p'
    if sort_col != 'mlog10p':
        df[sort_col] = pd.to_numeric(df[sort_col], errors='coerce')
    df['abs_sort'] = df[sort_col].abs()
    res = df.nlargest(top_n, 'abs_sort')
    res = res.sort_values(sort_col, ascending=True)
    return res

def load_ext_volcano_data(f1):
    df = load_ext_data_basic(f1)
    if df is None: return None, None
    ratio_col = 'OddsRatio' if 'OddsRatio' in df.columns else ('HazardRatio' if 'HazardRatio' in df.columns else 'FoldChange')
    if ratio_col not in df.columns: return None, None
    df[ratio_col] = pd.to_numeric(df[ratio_col], errors='coerce')
    df['log2Ratio'] = np.log2(df[ratio_col].replace(0, np.nan)) 
    return df, ratio_col

def load_knockdown_file(file_path):
    if not os.path.exists(file_path):
        return None, None
    try:
        try: df = pd.read_csv(file_path, index_col=0, encoding='utf-8-sig', low_memory=False)
        except: df = pd.read_csv(file_path, index_col=0, encoding='cp950', low_memory=False)

        if 'Subtype' not in df.index:
            print("Error: 'Subtype' row not found in knockdown file.")
            return None, None
        
        subtypes = df.loc['Subtype']
        
        gene_data = df.drop([idx for idx in ['CellLine', 'Subtype'] if idx in df.index])
        
        gene_data_transposed = gene_data.T
        
        gene_data_transposed = gene_data_transposed.apply(pd.to_numeric, errors='coerce')
        
        return gene_data_transposed, subtypes
    except Exception as e:
        print(f"Error reading knockdown file: {e}")
        return None, None

# ================= 2.1 Drug Target Annotation Helpers =================
def resolve_file_path(filename):
    if os.path.exists(filename):
        return filename
    try:
        alt = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
        if os.path.exists(alt):
            return alt
    except:
        pass
    return filename

def build_gene_drug_map(filepath="interactions.csv"):
    full_path = resolve_file_path(filepath)
    if not os.path.exists(full_path):
        print(f"\n[System] Drug target file not found: {full_path}")
        return {}
    try:
        # Use engine='python' and sep=None to auto-detect delimiters, avoiding issues with tabs or extra spaces
        try:
            df = pd.read_csv(full_path, encoding='utf-8-sig', sep=None, engine='python')
        except:
            try:
                df = pd.read_csv(full_path, encoding='cp950', sep=None, engine='python')
            except:
                with open(full_path, 'r', encoding='utf-8', errors='replace') as f:
                    df = pd.read_csv(f, sep=None, engine='python')
        
        # Ultimate cleanup: remove all non-alphanumeric characters from column names and convert to lowercase
        df.columns = df.columns.astype(str).str.replace(r'[^a-zA-Z0-9]', '', regex=True).str.lower()
        
        mapping = {}
        has_approved = 'approved' in df.columns
        for _, row in df.iterrows():
            gene = str(row.get('genename', '')).strip().upper()
            drug = str(row.get('drugname', '')).strip()
            if not gene or gene in ['NAN', 'NULL', 'NONE'] or not drug or drug in ['NAN', 'NULL', 'NONE']: continue
            drug_label = drug
            if has_approved:
                approved = str(row.get('approved', '')).strip().lower()
                if approved in ['true', '1', 'yes', 'y']:
                    drug_label += "(Approved)"
            if gene not in mapping: mapping[gene] = set()
            mapping[gene].add(drug_label)
        return {g: ", ".join(sorted(list(ds))) for g, ds in mapping.items()}
    except Exception as e:
        print(f"\n[Error] Exception while reading drug data: {e}")
        return {}

def build_pathway_gene_map(filepath="Full_Gene_List.csv"):
    full_path = resolve_file_path(filepath)
    if not os.path.exists(full_path):
        return {}
    try:
        try: df = pd.read_csv(full_path, encoding='utf-8-sig')
        except: df = pd.read_csv(full_path, encoding='cp950')
        mapping = {}
        for _, row in df.iterrows():
            pathway_raw = str(row.iloc[0]).strip().upper()
            genes_raw = str(row.iloc[2]).strip().upper()
            if not pathway_raw or pathway_raw == 'NAN' or pd.isna(genes_raw) or genes_raw == 'NAN': continue
            genes = [g.strip() for g in genes_raw.split(',') if g.strip()]
            mapping[pathway_raw] = genes
            pathway_clean = pathway_raw.split('=')[0].strip()
            mapping[pathway_clean] = genes
        return mapping
    except Exception as e:
        return {}

def annotate_top_tables(df, type_name, prefix, gene_to_drugs, pathway_to_genes):
    df = df.copy()
    col_name = f"{prefix}_Drug_Targets"
    if df.empty:
        df[col_name] = []
        return df
    
    name_col = df.columns[0] 
    if not gene_to_drugs:
        df[col_name] = "Target Data Not Loaded"
        return df
        
    if type_name == "Gene":
        df[col_name] = df[name_col].map(lambda x: gene_to_drugs.get(str(x).strip().upper(), "No Target"))
    else:
        def get_path_drugs(p_name):
            p_name_clean = str(p_name).strip().upper()
            genes = pathway_to_genes.get(p_name_clean, [])
            if not genes:
                for k, v in pathway_to_genes.items():
                    if p_name_clean.startswith(k) or k.startswith(p_name_clean):
                        genes = v; break
            drugs_found = [f"{g}: [{gene_to_drugs[g]}]" for g in genes if g in gene_to_drugs]
            return " | ".join(drugs_found) if drugs_found else "No Target"
        df[col_name] = df[name_col].map(get_path_drugs)
    
    df[col_name] = df[col_name].replace("", "No Target")
    return df

# ================= 3. Core Calculation Functions =================
def calculate_pearson(vec1, vec2):
    """
    Calculate Pearson correlation after correctly aligning samples by index.

    This is important when vec1 and vec2 come from different matrices,
    for example:
        Gene_expression.csv index order
        vs
        Signaling_pathway.csv index order

    Without index alignment, the same 30 cell lines may be paired in the wrong order.
    """

    s1 = pd.Series(vec1).copy()
    s2 = pd.Series(vec2).copy()

    s1 = pd.to_numeric(s1, errors='coerce')
    s2 = pd.to_numeric(s2, errors='coerce')

    # Critical fix:
    # align by sample / cell line ID before correlation
    df_pair = pd.concat(
        [s1.rename("x"), s2.rename("y")],
        axis=1,
        join="inner"
    )

    df_pair = df_pair.dropna()

    if df_pair.shape[0] < 3:
        return None, None, 0

    c1 = df_pair["x"]
    c2 = df_pair["y"]

    if np.std(c1) == 0 or np.std(c2) == 0:
        return None, None, len(c1)

    try:
        r, p = stats.pearsonr(c1, c2)
        return r, p, len(c1)
    except:
        return None, None, len(c1)

def run_anova(data_df, group_col, value_col, hue_col=None):
    cols_to_keep = [group_col, value_col]
    if hue_col:
        cols_to_keep.append(hue_col)
        
    clean_df = data_df[cols_to_keep].copy()
    clean_df[value_col] = pd.to_numeric(clean_df[value_col], errors='coerce')
    clean_df = clean_df.dropna()
    
    groups = []
    group_stats = []
    
    if hue_col:
        combined_col = f"{group_col}_{hue_col}"
        clean_df[combined_col] = clean_df[group_col].astype(str) + " & " + clean_df[hue_col].astype(str)
        unique_cats = clean_df[combined_col].unique()
        target_col = combined_col
    else:
        unique_cats = clean_df[group_col].unique()
        target_col = group_col
        
    for cat in unique_cats:
        vals = clean_df[clean_df[target_col] == cat][value_col]
        groups.append(vals)
        group_stats.append({'Group': cat, 'N': len(vals), 'Mean': vals.mean(), 'Std': vals.std()})
        
    if len(groups) < 2: 
        return None, pd.DataFrame(group_stats)
        
    try:
        if len(groups) == 2:
            t_stat, p_val = stats.ttest_ind(groups[0], groups[1], equal_var=False)
            return p_val, pd.DataFrame(group_stats)
        else:
            f_stat, p_val = stats.f_oneway(*groups)
            return p_val, pd.DataFrame(group_stats)
    except: 
        return None, pd.DataFrame(group_stats)

def calculate_signature_score(gene_list, df_gene, method='1'):
    valid_genes = [g for g in gene_list if g in df_gene.columns]
    if not valid_genes: return None, []
    sub_df = df_gene[valid_genes].copy()
    if method == '3':
        rank_df = sub_df.rank(pct=True)
        score = rank_df.mean(axis=1)
    else:
        z_df = (sub_df - sub_df.mean()) / (sub_df.std() + 1e-9)
        score = z_df.mean(axis=1)
        if method == '2': score = score.clip(lower=-3, upper=3)
    return score, valid_genes

def get_data_by_name(name, df_gene, df_pathway):
    name = name.strip()
    if name in df_gene.columns: return df_gene[name], "Gene"
    elif name in df_pathway.columns: return df_pathway[name], "Pathway"
    if '/' in name:
        parts = name.split('/')
        if len(parts) == 2:
            n1, n2 = parts[0].strip(), parts[1].strip()
            d1, _ = get_data_by_name(n1, df_gene, df_pathway)
            d2, _ = get_data_by_name(n2, df_gene, df_pathway)
            if d1 is not None and d2 is not None:
                return d1 / d2.replace(0, 1e-9), "Virtual_Ratio"
    if ' - ' in name:
        parts = name.split(' - ')
        if len(parts) == 2:
            n1, n2 = parts[0].strip(), parts[1].strip()
            d1, _ = get_data_by_name(n1, df_gene, df_pathway)
            d2, _ = get_data_by_name(n2, df_gene, df_pathway)
            if d1 is not None and d2 is not None:
                return d1 - d2, "Virtual_Diff"
    return None, None

def detect_survival_pairs(df_clinical):
    """
    Automatically detect possible survival time and status column pairs in clinical data.
    For example: OS_Months and OS_Status, DFS_Time and DFS_Event.
    """
    survival_pairs = []
    
    # 收集所有可能的 time 和 status 欄位
    all_time_cols = [col for col in df_clinical.columns if 'months' in col.lower() or 'time' in col.lower()]
    all_status_cols = [col for col in df_clinical.columns if 'status' in col.lower() or 'event' in col.lower()]

    # Extract prefixes and try to pair them
    prefixes = set()
    for col in all_time_cols + all_status_cols:
        match = re.match(r'([a-zA-Z]+)[_ ]*(months|time|status|event)', col, re.IGNORECASE)
        if match:
            prefixes.add(match.group(1))

    for prefix in sorted(list(prefixes)):
        time_col_found = next((col for col in all_time_cols if col.lower().startswith(prefix.lower()) and ('months' in col.lower() or 'time' in col.lower())), None)
        status_col_found = next((col for col in all_status_cols if col.lower().startswith(prefix.lower()) and ('status' in col.lower() or 'event' in col.lower())), None)

        if time_col_found and status_col_found:
            survival_pairs.append((time_col_found, status_col_found))
            
    # If no pairs are detected, but default OS_Months and OS_Status exist, add them
    if not survival_pairs and 'OS_Months' in df_clinical.columns and 'OS_Status' in df_clinical.columns:
        survival_pairs.append(('OS_Months', 'OS_Status'))
    
    # Ensure uniqueness and sort
    survival_pairs = sorted(list(set(survival_pairs)))
    return survival_pairs

# ================= 4. Plotting and Helper Functions =================
def clean_filename(name):
    return re.sub(r'[\\/*?:"<>|]', "_", str(name))

def get_desktop_path():
    base_desktop = os.path.join(os.path.expanduser("~"), "Desktop")
    out_folder = os.path.join(base_desktop, SESSION_FOLDER_NAME)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    return out_folder

def setup_plot_style():
    plt.close('all') 
    plt.rcdefaults() 
    sns.set(style="ticks") 
    warnings.filterwarnings("ignore", message=".*Glyph.*")
    plt.rcParams.update({
        'font.sans-serif': ['Microsoft JhengHei', 'SimHei', 'Arial', 'sans-serif'],
        'axes.unicode_minus': False,
        'figure.dpi': IMG_DPI,
        'axes.titlesize': 20,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'legend.title_fontsize': 16,
        'font.size': 14
    })

def is_clinical_noise(s):
    return bool(re.search(r'(^|[_ \-\.])(x|unknown|unk|na|n/a|not reported|not_reported|missing)($|[_ \-\.])', str(s).lower()))

def bio_sort_key(s):
    s_str = str(s).lower()
    if is_clinical_noise(s):
        prio = 4 
    elif re.search(r'(^|[_ \-\.])(pre|baseline|pretreatment|wt|wildtype|wild type|wild-type|normal|control|negative)($|[_ \-\.])', s_str):
        prio = 0
    elif re.search(r'(^|[_ \-\.])(post|posttreatment|mut|mutant|mutation|tumor|cancer|positive)($|[_ \-\.])', s_str):
        prio = 2
    else:
        prio = 1
        
    alphanum = [f"{float(t):010.4f}" if t.replace('.','',1).isdigit() else t for t in re.split('([0-9.]+)', s_str)]
    return (prio, alphanum)

def run_patient_waterfall(target_name, target_data, df_clinical, subtype_col, subtypes_to_run=None, palette='Set1'):
    print(f"\n>>> Plotting patient expression waterfall for {target_name}...")
    if subtypes_to_run is not None and df_clinical is not None and subtype_col:
        if None in subtypes_to_run: valid_idx = df_clinical.index
        else: valid_idx = df_clinical[df_clinical[subtype_col].isin(subtypes_to_run)].index
        data = target_data.loc[target_data.index.intersection(valid_idx)].dropna()
    else:
        data = target_data.dropna()

    if len(data) < 2: 
        print("Insufficient data to plot.")
        return
        
    centered_data = data - data.mean()
    sorted_data = centered_data.sort_values(ascending=True)
    
    setup_plot_style()
    plt.figure(figsize=(max(10, len(sorted_data)*0.02), 7))
    x_pos = np.arange(len(sorted_data))
    
    if subtype_col and df_clinical is not None and subtype_col in df_clinical.columns:
        subtypes = df_clinical.loc[sorted_data.index, subtype_col].fillna("Unknown")
        unique_subtypes = list(subtypes.unique())
        
        def st_sort_key(s):
            return 1 if is_clinical_noise(s) else 0
        unique_subtypes.sort(key=st_sort_key)

        palette_colors = sns.color_palette(palette, n_colors=len(unique_subtypes))
        color_map = {}
        p_idx = 0
        for st in unique_subtypes:
            if is_clinical_noise(st):
                color_map[st] = "#B0B0B0"
            else:
                color_map[st] = palette_colors[p_idx]
                p_idx += 1
                
        colors = [color_map[st] for st in subtypes]
        
        plt.bar(x_pos, sorted_data, color=colors, width=1.0, edgecolor='none')
        import matplotlib.patches as mpatches
        patches = [mpatches.Patch(color=color_map[st], label=st) for st in unique_subtypes]
        plt.legend(handles=patches, title=subtype_col, loc='upper right', title_fontsize='18', fontsize='16')
    else:
        colors = ['#d62728' if v >= 0 else '#1f77b4' for v in sorted_data]
        plt.bar(x_pos, sorted_data, color=colors, width=1.0, edgecolor='none')
        
    plt.axhline(0, color='black', linewidth=1)
    plt.title(f"Patient Waterfall: {target_name}\n(Mean-centered)", fontsize=22, fontweight='bold')
    plt.xlabel(f"Patients (N={len(sorted_data)})", fontsize=18)
    plt.ylabel("Relative Value (Mean-centered)", fontsize=18)
    plt.xticks([]) 
    plt.tick_params(axis='y', labelsize=16)
    plt.tight_layout()
    fname = f"Waterfall_Patient_{clean_filename(target_name)}.{IMG_EXT}"
    plt.savefig(os.path.join(get_desktop_path(), fname), dpi=IMG_DPI)
    plt.close()
    print(f"-> [Done] Patient waterfall plot saved: {fname}")

def show_scatter_plot(vec1, vec2, name1, name2, r_val, p_val, n_val):
    setup_plot_style()
    df_plot = pd.DataFrame({name1: vec1, name2: vec2}).dropna()
    fig, ax = plt.subplots(figsize=(8, 7))
    _draw_scatter(ax, df_plot, name1, name2, r_val, p_val, n_val, "Correlation Analysis")
    fname = f"Quick_Plot_{clean_filename(name1)}_vs_{clean_filename(name2)}.{IMG_EXT}"
    try: plt.savefig(os.path.join(get_desktop_path(), fname), dpi=IMG_DPI)
    except: pass
    try: plt.show(); plt.close()
    except: pass

def save_subtype_scatter(df_data, name1, name2, output_folder, subtype_name, r_val, p_val, n_val):
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 7))
    title = f"{name1} vs {name2}\n({subtype_name})"
    _draw_scatter(ax, df_data, name1, name2, r_val, p_val, n_val, title)
    fname = f"Corr_{subtype_name}_{clean_filename(name1)}_vs_{clean_filename(name2)}.{IMG_EXT}"
    plt.savefig(os.path.join(output_folder, fname), dpi=IMG_DPI)
    plt.close()

# [v81.1 Fix] Matplotlib requires plural arguments: linewidths & edgecolors
def _draw_scatter(ax, df, x, y, r, p, n, title):
    df = df.copy()
    df[x] = pd.to_numeric(df[x], errors='coerce')
    df[y] = pd.to_numeric(df[y], errors='coerce')
    df = df.dropna()
    
    x_nunique = df[x].nunique()
    y_nunique = df[y].nunique()
    x_range = df[x].max() - df[x].min() if x_nunique > 1 else 1
    y_range = df[y].max() - df[y].min() if y_nunique > 1 else 1
    
    x_jit = x_range * 0.03 if x_nunique < 15 else None
    y_jit = y_range * 0.03 if y_nunique < 15 else None

    alpha_val = 0.3 if (x_jit or y_jit) else 0.5
    s_val = 30 if (x_jit or y_jit) else 45

    sns.regplot(data=df, x=x, y=y, 
                x_jitter=x_jit, y_jitter=y_jit,
                scatter_kws={'s': s_val, 'alpha': alpha_val, 'edgecolors': 'white', 'linewidths': 0.5}, 
                line_kws={'color': '#d62728', 'linewidth': 2.5}, ax=ax)
                
    ax.set_title(title, fontsize=22, fontweight='bold', pad=40)
    ax.set_xlabel(x, fontsize=20)
    ax.set_ylabel(y, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)

    p_color = '#d62728' if p < 0.05 else 'black'
    p_weight = 'bold' if p < 0.05 else 'normal'
    r_text = f"$R$ = {r:.3f}, $N$ = {n}"
    p_text = f"p {format_pval(p, with_stars=True)}"

    ax.text(0.97, 0.97, r_text, transform=ax.transAxes, ha='right', va='top', fontsize=18, bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.7))
    ax.text(0.5, 1.02, p_text, transform=ax.transAxes, ha='center', va='bottom', fontsize=18, color=p_color, fontweight=p_weight)

def save_clinical_boxplot(df_merged, feature, target_name, output_folder, subtype_name, p_val, hue_col=None, palette='Set1'):
    setup_plot_style()
    cols_to_keep = [feature, 'Target_Value']
    if hue_col: cols_to_keep.append(hue_col)
    plot_data = df_merged[cols_to_keep].copy()
    plot_data['Target_Value'] = pd.to_numeric(plot_data['Target_Value'], errors='coerce')
    plot_data = plot_data.dropna()
    if plot_data.empty or len(plot_data[feature].unique()) < 1: return

    # Force to string to prevent Seaborn palette key errors with float/int
    plot_data[feature] = plot_data[feature].astype(str)
    if hue_col:
        plot_data[hue_col] = plot_data[hue_col].astype(str)

    try:
        unique_vals = list(plot_data[feature].unique())
        
        sorted_cats = sorted(unique_vals, key=bio_sort_key)
        plot_data[feature] = pd.Categorical(plot_data[feature], categories=sorted_cats, ordered=True)
        
        sorted_hues = None
        dynamic_palette = {}
        if hue_col:
            hue_vals = list(plot_data[hue_col].unique())
            sorted_hues = sorted(hue_vals, key=bio_sort_key)
            plot_data[hue_col] = pd.Categorical(plot_data[hue_col], categories=sorted_hues, ordered=True)
            n_hues = len(sorted_hues)
            colors = sns.color_palette(palette, n_colors=n_hues).as_hex()
            dynamic_palette = {c: colors[i] for i, c in enumerate(sorted_hues)}
        else:
            valid_cats = [c for c in sorted_cats if not is_clinical_noise(c)]
            noise_cats = [c for c in sorted_cats if is_clinical_noise(c)]
            
            n_valid = len(valid_cats)
            if n_valid == 1:
                dynamic_palette[valid_cats[0]] = "#1f77b4" 
            elif n_valid == 2:
                dynamic_palette[valid_cats[0]] = "#1f77b4" 
                dynamic_palette[valid_cats[1]] = "#d62728" 
            elif n_valid > 2:
                colors = sns.color_palette("coolwarm", n_colors=n_valid).as_hex()
                for i, c in enumerate(valid_cats):
                    dynamic_palette[c] = colors[i]
            for c in noise_cats:
                dynamic_palette[c] = "#B0B0B0"
            
    except Exception as e:
        dynamic_palette = palette

    n_cats_total = len(plot_data[feature].unique())
    fig_w = max(8, n_cats_total * (0.8 if hue_col else 0.5))
    fig, ax = plt.subplots(figsize=(fig_w, 7))

    if hue_col:
        test_type = 'Welch t-test' if len(plot_data[feature].unique()) * len(plot_data[hue_col].unique()) == 2 else 'ANOVA'
    else:
        test_type = 'Welch t-test' if len(plot_data[feature].unique()) == 2 else 'ANOVA'
        
    title_suffix = f"by {feature}" + (f" & {hue_col}" if hue_col else "") + f"\n({subtype_name})"

    sns.violinplot(data=plot_data, x=feature, y='Target_Value', hue=hue_col if hue_col else None,
                   order=sorted_cats, hue_order=sorted_hues if hue_col else None,
                   palette=dynamic_palette, inner=None, linewidth=1, alpha=0.4, ax=ax)

    sns.boxplot(data=plot_data, x=feature, y='Target_Value', hue=hue_col if hue_col else None,
                order=sorted_cats, hue_order=sorted_hues if hue_col else None,
                width=0.4 if hue_col else 0.15, color='#333333',
                showfliers=False, showmeans=False, dodge=True, 
                boxprops={'facecolor':'none', 'edgecolor':'#333333', 'linewidth':2}, 
                medianprops={'color':'#FFFF33', 'linewidth':3}, 
                whiskerprops={'color':'#333333', 'linewidth':1.5}, capprops={'color':'#333333', 'linewidth':1.5}, ax=ax)
        
    sns.stripplot(data=plot_data, x=feature, y='Target_Value', hue=hue_col if hue_col else None,
                  order=sorted_cats, hue_order=sorted_hues if hue_col else None,
                  color='#333333', size=3, alpha=0.4, jitter=True, dodge=True, ax=ax)

    if hue_col:
        handles, labels = ax.get_legend_handles_labels()
        n_h = len(sorted_hues)
        ax.legend(handles[:n_h], labels[:n_h], title=hue_col, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

    ax.set_title(f"{target_name} {title_suffix}", fontsize=22, fontweight='bold', pad=40)
    p_color = '#d62728' if p_val is not None and p_val < 0.05 else 'black'
    p_text = f"Overall {test_type} p {format_pval(p_val, with_stars=True)}"
    ax.text(0.5, 1.02, p_text, transform=ax.transAxes, ha='center', va='bottom', fontsize=18, fontweight='bold', color=p_color)

    ax.set_xlabel(feature, fontsize=16); ax.set_ylabel("Expression / Score", fontsize=16)
    ax.tick_params(axis='x', rotation=45)
    fig.tight_layout()
    fname_suffix = f"_{clean_filename(hue_col)}" if hue_col else ""
    fig.savefig(os.path.join(output_folder, f"{subtype_name}_{clean_filename(feature)}{fname_suffix}.{IMG_EXT}"), dpi=IMG_DPI)
    plt.close(fig)

def save_high_low_boxplot(df_plot, group_name, compare_name, p_val, output_folder, subtype_name):
    setup_plot_style() # This resets global settings, might be better inside the calling function
    fig, ax = plt.subplots(figsize=(7, 7))
    order = ["Low", "High"]
    locked_palette = {"Low": "#1f77b4", "High": "#d62728"}
    
    title_prefix = f"{compare_name}\nby {group_name} Status ({subtype_name})"

    sns.violinplot(data=df_plot, x='Group', y='Value', order=order,
                   palette=locked_palette, inner=None, linewidth=1, alpha=0.4, ax=ax)

    sns.boxplot(data=df_plot, x='Group', y='Value', order=order,
                width=0.15, color='#333333',
                showfliers=False, showmeans=False, dodge=False, 
                boxprops={'facecolor':'none', 'edgecolor':'#333333', 'linewidth':2},
                medianprops={'color':'#FFFF33', 'linewidth':3},
                whiskerprops={'color':'#333333', 'linewidth':1.5}, capprops={'color':'#333333', 'linewidth':1.5}, ax=ax)
    
    sns.stripplot(data=df_plot, x='Group', y='Value', order=order, color='#333333', size=3, alpha=0.4, jitter=True, ax=ax)
    
    ax.set_title(title_prefix, fontsize=22, fontweight='bold', pad=40)
    p_color = '#d62728' if p_val < 0.05 else 'black'
    p_text = f"Welch t-test p {format_pval(p_val, with_stars=True)}"
    ax.text(0.5, 1.02, p_text, transform=ax.transAxes, ha='center', va='bottom', fontsize=18, fontweight='bold', color=p_color)

    ax.set_xlabel(f"{group_name} Status", fontsize=20)
    ax.set_ylabel("Expression / Ratio", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    fig.tight_layout()
    
    fname = f"HighLow_{clean_filename(group_name)}_vs_{clean_filename(compare_name)}_{clean_filename(subtype_name)}.{IMG_EXT}"
    fig.savefig(os.path.join(output_folder, fname), dpi=IMG_DPI)
    plt.close(fig)
def run_high_low_comparison(group_name, group_data, compare_dict, df_clinical, subtype_col, subtypes_to_run=None):
    print(f"\n>>> Running high/low group comparison (grouped by median of {group_name})...")
    
    img_folder = os.path.join(get_desktop_path(), f"HighLow_Plots_{clean_filename(group_name)}")
    if not os.path.exists(img_folder): os.makedirs(img_folder)
    
    groups = get_analysis_indices(df_clinical, group_data.index, subtype_col, subtypes_to_run)
            
    results = []
    for grp_name, grp_index in groups:
        valid_idx = grp_index.intersection(group_data.dropna().index)
        if len(valid_idx) < 2: continue
        
        sub_group_data = group_data.loc[valid_idx]
        median_val = sub_group_data.median()
        group_labels = sub_group_data.apply(lambda x: "High" if x >= median_val else "Low")
        
        for comp_name, comp_data in compare_dict.items():
            comp_idx = valid_idx.intersection(comp_data.dropna().index)
            if len(comp_idx) < 2: continue
            
            df_plot = pd.DataFrame({
                'Group': group_labels.loc[comp_idx],
                'Value': comp_data.loc[comp_idx]
            })
            
            high_vals = df_plot[df_plot['Group'] == 'High']['Value']
            low_vals = df_plot[df_plot['Group'] == 'Low']['Value']
            if len(high_vals) < 1 or len(low_vals) < 1: continue
            
            if len(high_vals) >= 2 and len(low_vals) >= 2:
                t_stat, p_val = stats.ttest_ind(high_vals, low_vals, equal_var=False)
            else:
                t_stat, p_val = np.nan, np.nan
            save_high_low_boxplot(df_plot, group_name, comp_name, p_val, img_folder, grp_name)
            
            results.append({
                'Scope': grp_name,
                'Group_By': group_name,
                'Compare_Target': comp_name,
                'N_High': len(high_vals),
                'Mean_High': high_vals.mean(),
                'N_Low': len(low_vals),
                'Mean_Low': low_vals.mean(),
                'T_Stat': t_stat,
                'P_Value': p_val
            })
            
            try:
                raw_fname = f"RawData_HighLow_{clean_filename(group_name)}_vs_{clean_filename(comp_name)}_{clean_filename(grp_name)}.xlsx"
                df_plot.to_excel(os.path.join(img_folder, raw_fname), index_label="Patient_ID")
            except: pass

    if results:
        df_res = pd.DataFrame(results)
        res_fname = f"HighLow_Summary_{clean_filename(group_name)}.xlsx"
        
        cols_to_show = ['Scope', 'Group_By', 'Compare_Target', 'N_High', 'Mean_High', 'N_Low', 'Mean_Low', 'T_Stat', 'P_Value']
        if FDR_ENABLED:
            # [FDR] Apply FDR correction to P_Values within each Scope
            df_res['Adjusted_P_Value'] = df_res.groupby('Scope')['P_Value'].transform(lambda p: multipletests(p.dropna(), method='fdr_bh')[1] if not p.dropna().empty else np.nan)
            cols_to_show.append('Adjusted_P_Value')
        df_res = df_res[cols_to_show]

        df_res.to_excel(os.path.join(img_folder, res_fname), index=False)
        
def run_clinical_analysis(target_name, target_data, df_clinical, subtype_col, plot_enabled, selected_features=None, hue_col=None, subtypes_to_run=None, palette='Set1'):
    print(f"\n>>> Analyzing clinical features for {target_name} (Box+Violin)...")
    merged_df = df_clinical.copy()
    merged_df['Target_Value'] = target_data
    try:
        raw_fname = f"RawData_Clinical_{clean_filename(target_name)}.xlsx"
        merged_df.to_excel(os.path.join(get_desktop_path(), raw_fname), index_label="Patient_ID")
    except: pass
    img_folder = os.path.join(get_desktop_path(), f"Clinical_Plots_{clean_filename(target_name)}")
    if plot_enabled and not os.path.exists(img_folder): os.makedirs(img_folder)
    
    results, cols = [], []
    if selected_features:
        for c in selected_features:
            if c in df_clinical.columns and c != 'Target_Value':
                cols.append(c)
        print(f"  -> [System] User specified features for analysis: {', '.join(cols)}")
    else:
        for c in df_clinical.columns:
            if c == subtype_col or c == 'Target_Value': continue
            c_lower = str(c).lower()
            # [Precise Debug] Exclude only true ID columns, but allow 'aneuploIDy'
            if c_lower in ['id', 'patient_id', 'sample_id', 'patientid', 'sampleid'] or c_lower.endswith('_id') or c_lower.startswith('id_'): continue
            if any(k in c_lower for k in ['time', 'month', 'day', 'survival', 'os_', 'rfs_', 'dfs_', 'pfs_']): continue
            if pd.api.types.is_numeric_dtype(df_clinical[c]) and len(df_clinical[c].dropna().unique()) > 150: continue
            cols.append(c)
        print(f"  -> [System] Automatically detected features for analysis: {', '.join(cols)}")

    if not cols:
        print("  -> [Warning] No analyzable clinical features found!")
        return

    groups = []
    if subtypes_to_run is None or (len(subtypes_to_run) == 1 and subtypes_to_run[0] is None):
        groups.append(('All_Patients', merged_df))
        if subtype_col and subtype_col in merged_df.columns:
            groups += [(f"Sub_{st}", merged_df[merged_df[subtype_col]==st]) for st in merged_df[subtype_col].dropna().unique()]
    else:
        for st in subtypes_to_run:
            if st is None: groups.append(('All_Patients', merged_df))
            elif subtype_col and subtype_col in merged_df.columns: groups.append((f"Sub_{st}", merged_df[merged_df[subtype_col]==st]))
        
    for grp_name, grp_df in groups:
        for feat in cols:
            if feat == hue_col: continue # 防止自己跟自己比較
            p_val, stats_df = run_anova(grp_df, feat, 'Target_Value', hue_col=hue_col)
            if plot_enabled: save_clinical_boxplot(grp_df, feat, target_name, img_folder, grp_name, p_val, hue_col=hue_col, palette=palette)
            for _, row in stats_df.iterrows(): results.append({'Scope':grp_name, 'Factor':feat, 'Hue_Factor': hue_col, 'Group':row['Group'], 'N':row['N'], 'Mean':row['Mean'], 'Std':row['Std'], 'P':p_val})
    
    if results:
        df_res = pd.DataFrame(results)
        if FDR_ENABLED:
            # [FDR] Apply FDR correction to P-values within each Scope
            df_res['Adjusted_P_Value'] = df_res.groupby('Scope')['P'].transform(lambda p: multipletests(p.dropna().unique(), method='fdr_bh')[1][pd.factorize(p.dropna())[0]] if not p.dropna().empty else np.nan)
        df_res.to_excel(os.path.join(get_desktop_path(), f"Clinical_{clean_filename(target_name)}.xlsx"), index=False)

def run_survival_analysis(target_name, target_data, df_clinical, subtype_col, time_col, status_col, subtypes_to_run=None):
    print(f"\n>>> Analyzing survival curve for {target_name} (KM Plot)...")
    merged_df = df_clinical.copy()
    merged_df['Target_Value'] = target_data
    if time_col not in merged_df.columns or status_col not in merged_df.columns:
        print(f"Error: Survival columns not found."); return
    
    img_folder = os.path.join(get_desktop_path(), f"Survival_{clean_filename(target_name)}")
    if not os.path.exists(img_folder): os.makedirs(img_folder)
    results = []
    groups = get_analysis_groups(merged_df, subtype_col, subtypes_to_run)
    
    for grp_name, grp_df in groups:
        sub_df = grp_df[['Target_Value', time_col, status_col]].copy()
        sub_df['Target_Value'] = pd.to_numeric(sub_df['Target_Value'], errors='coerce')
        sub_df = sub_df.dropna()
        if len(sub_df) < 2: continue
        
        median_val = sub_df['Target_Value'].median()
        high_group = sub_df[sub_df['Target_Value'] >= median_val]
        low_group = sub_df[sub_df['Target_Value'] < median_val]
        if len(high_group) < 1 or len(low_group) < 1: continue
        
        try:
            results_lr = logrank_test(high_group[time_col], low_group[time_col], event_observed_A=high_group[status_col], event_observed_B=low_group[status_col])
            p_val = results_lr.p_value
        except:
            p_val = np.nan
        
        cph_df = sub_df.copy()
        cph_df['Group'] = (cph_df['Target_Value'] >= median_val).astype(int) 
        try:
            cph = CoxPHFitter(penalizer=0.05)
            cph.fit(cph_df[['Group', time_col, status_col]], duration_col=time_col, event_col=status_col)
            hr = cph.hazard_ratios_['Group']
            ci_lower = cph.confidence_intervals_.iloc[0, 0]
            ci_upper = cph.confidence_intervals_.iloc[0, 1]
            if hr > 100 or hr < 0.01:
                hr_str = "HR = N/A (Separation)"
            else:
                hr_str = f"HR = {hr:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f})"
        except:
            hr_str = "HR = N/A"

        fig, ax = plt.subplots(figsize=(8, 7))
        kmf = KaplanMeierFitter()
        kmf.fit(high_group[time_col], event_observed=high_group[status_col], label=f"High (N={len(high_group)})")
        kmf.plot_survival_function(ax=ax, color="#d62728", ci_show=False)
        kmf.fit(low_group[time_col], event_observed=low_group[status_col], label=f"Low (N={len(low_group)})")
        kmf.plot_survival_function(ax=ax, color="#1f77b4", ci_show=False)
        
        ax.set_title(f"{target_name} Survival ({grp_name})\n{hr_str}", fontsize=22, fontweight='bold', pad=40)
        p_color = '#d62728' if p_val < 0.05 else 'black'
        p_text = f"Log-rank p {format_pval(p_val, with_stars=True)}"
        ax.text(0.5, 1.02, p_text, transform=ax.transAxes, ha='center', va='bottom', fontsize=18, fontweight='bold', color=p_color)
        ax.set_xlabel(f"Time ({time_col})", fontsize=20); ax.set_ylabel("Probability", fontsize=20); ax.set_ylim(0, 1.05)
        ax.tick_params(axis='both', which='major', labelsize=16); ax.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(os.path.join(img_folder, f"KM_{grp_name}_{clean_filename(target_name)}_{clean_filename(time_col)}_{clean_filename(status_col)}.{IMG_EXT}"), dpi=IMG_DPI)
        plt.close(fig)
        
        results.append({'Scope': grp_name, 'Median': median_val, 'N_High': len(high_group), 'N_Low': len(low_group), 'LogRank_P': p_val, 'Time_Col': time_col, 'Status_Col': status_col})
    
    if results:
        pd.DataFrame(results).to_excel(os.path.join(get_desktop_path(), f"Survival_Stats_{clean_filename(target_name)}_{clean_filename(time_col)}_{clean_filename(status_col)}.xlsx"), index=False)

def run_dual_survival_analysis(name1, data1, name2, data2, df_clinical, subtype_col, time_col, status_col, subtypes_to_run=None):
    print(f"  -> Plotting dual-target survival curve for {name1} + {name2}...")
    merged_df = df_clinical.copy()
    merged_df['T1'] = pd.to_numeric(data1, errors='coerce')
    merged_df['T2'] = pd.to_numeric(data2, errors='coerce')
    
    if time_col not in merged_df.columns or status_col not in merged_df.columns:
        print(f"     [Warning] Survival columns {time_col} or {status_col} not found.")
        return

    img_folder = os.path.join(get_desktop_path(), f"Dual_Survival_{clean_filename(name1)}")
    if not os.path.exists(img_folder): os.makedirs(img_folder)

    groups = get_analysis_groups(merged_df, subtype_col, subtypes_to_run)

    results = []
    for grp_name, grp_df in groups:
        sub_df = grp_df[['T1', 'T2', time_col, status_col]].dropna()
        if len(sub_df) < 2: continue

        med1 = sub_df['T1'].median()
        med2 = sub_df['T2'].median()

        def assign_group(row):
            g1 = "High" if row['T1'] >= med1 else "Low"
            g2 = "High" if row['T2'] >= med2 else "Low"
            return f"{g1}_{g2}"

        sub_df['Combined'] = sub_df.apply(assign_group, axis=1)
        if sub_df['Combined'].nunique() < 2: continue

        # 1. Calculate Overall P-value
        try:
            res_lr = multivariate_logrank_test(sub_df[time_col], sub_df['Combined'], sub_df[status_col])
            overall_p = res_lr.p_value
        except:
            overall_p = np.nan

        base_grp = "Low_Low"
        # 2. Dummy Variable Cox HR calculation for all samples (Baseline: Low_Low)
        hr_results = {}
        if sub_df['Combined'].nunique() > 1 and base_grp in sub_df['Combined'].values:
            cph_df = sub_df[[time_col, status_col, 'Combined']].copy()
            cph_df = pd.get_dummies(cph_df, columns=['Combined'], drop_first=False)
            cph_df = cph_df.drop(columns=[f'Combined_{base_grp}'], errors='ignore')
            
            for col in cph_df.columns:
                cph_df[col] = cph_df[col].astype(float)
                
            try:
                cph = CoxPHFitter(penalizer=0.05)
                cph.fit(cph_df, duration_col=time_col, event_col=status_col)
                
                for grp in ["High_High", "High_Low", "Low_High"]:
                    col_name = f'Combined_{grp}'
                    if col_name in cph.summary.index:
                        hr = cph.hazard_ratios_[col_name]
                        ci_l = cph.confidence_intervals_.loc[col_name].iloc[0]
                        ci_u = cph.confidence_intervals_.loc[col_name].iloc[1]
                        p_cox = cph.summary.loc[col_name, 'p']
                        
                        if 0.01 <= hr <= 100:
                            hr_results[grp] = (hr, ci_l, ci_u, p_cox)
            except Exception:
                pass

        fig, ax = plt.subplots(figsize=(9, 7))
        kmf = KaplanMeierFitter()

        color_map = {"High_High": "#D62728", "Low_Low": "#1F77B4", "High_Low": "#FF7F0E", "Low_High": "#2CA02C"}
        label_base = {
            "High_High": f"{name1} High + {name2} High",
            "Low_Low": f"{name1} Low + {name2} Low",
            "High_Low": f"{name1} High + {name2} Low",
            "Low_High": f"{name1} Low + {name2} High"
        }

        for grp in ["High_High", "High_Low", "Low_High", "Low_Low"]:
            mask = sub_df['Combined'] == grp
            n_grp = mask.sum()
            if n_grp > 0:
                if grp == "Low_Low":
                    final_label = f"{label_base[grp]} (N={n_grp}, Ref)"
                else:
                    if grp in hr_results:
                        hr, ci_l, ci_u, p_cox = hr_results[grp]
                        stars = get_p_stars(p_cox)
                        sig_star = "" if stars == "n.s." else " " + stars
                        final_label = f"{label_base[grp]} (N={n_grp}, HR={hr:.2f}{sig_star})"
                    else:
                        final_label = f"{label_base[grp]} (N={n_grp})"
                        
                kmf.fit(sub_df.loc[mask, time_col], event_observed=sub_df.loc[mask, status_col], label=final_label)
                kmf.plot_survival_function(ax=ax, color=color_map[grp], ci_show=False, linewidth=2)

        ax.set_title(f"Dual-Target Survival: {name1} & {name2}\n({grp_name})", fontsize=22, fontweight='bold', pad=40)
        p_color = '#d62728' if overall_p < 0.05 else 'black'
        p_text = f"Overall p {format_pval(overall_p, with_stars=True)}"
        ax.text(0.5, 1.02, p_text, transform=ax.transAxes, ha='center', va='bottom', fontsize=18, fontweight='bold', color=p_color)
        ax.set_xlabel(f"Time ({time_col})", fontsize=20)
        ax.set_ylabel("Probability", fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_ylim(0, 1.05)
        ax.legend(loc='lower left', fontsize=13, framealpha=0.9, edgecolor='#CCCCCC')
        fig.tight_layout()
        fig.savefig(os.path.join(img_folder, f"DualKM_{grp_name}_{clean_filename(name1)}_vs_{clean_filename(name2)}_{clean_filename(time_col)}_{clean_filename(status_col)}.{IMG_EXT}"), dpi=IMG_DPI)
        plt.close(fig)

        res_item = {
            'Scope': grp_name,
            'Target1': name1, 'Target2': name2,
            'N_Total': len(sub_df),
            'N_HH': (sub_df['Combined'] == 'High_High').sum(),
            'N_HL': (sub_df['Combined'] == 'High_Low').sum(),
            'N_LH': (sub_df['Combined'] == 'Low_High').sum(),
            'N_LL': (sub_df['Combined'] == 'Low_Low').sum(),
            'Overall_P_Value': overall_p,
            'Time_Col': time_col, 'Status_Col': status_col
        }
        
        map_prefix = {"High_High": "HH", "High_Low": "HL", "Low_High": "LH"}
        for grp in ["High_High", "High_Low", "Low_High"]:
            prefix = map_prefix[grp]
            if grp in hr_results:
                hr, ci_l, ci_u, p_cox = hr_results[grp]
                res_item[f'HR ({prefix} vs LL)'] = hr
                res_item[f'95% CI ({prefix} vs LL)'] = f"{ci_l:.2f}-{ci_u:.2f}"
                res_item[f'P_Value ({prefix} vs LL)'] = p_cox
            else:
                res_item[f'HR ({prefix} vs LL)'] = np.nan
                res_item[f'95% CI ({prefix} vs LL)'] = ""
                res_item[f'P_Value ({prefix} vs LL)'] = np.nan
                
        results.append(res_item)

    if results:
        pd.DataFrame(results).to_excel(os.path.join(img_folder, f"DualKM_Stats_{clean_filename(name1)}_vs_{clean_filename(name2)}_{clean_filename(time_col)}_{clean_filename(status_col)}.xlsx"), index=False)

def run_subtype_correlation(name1, data1, name2, data2, df_clinical, subtype_col, plot_enabled, subtypes_to_run=None):
    print(f"\nAnalyzing {name1} vs {name2} ...")
    df_m = pd.DataFrame({name1:data1, name2:data2})
    if subtype_col: df_m = df_m.join(df_clinical[[subtype_col]], how='inner')
    
    try:
        raw_fname = f"RawData_SubCorr_{clean_filename(name1)}_vs_{clean_filename(name2)}.xlsx"
        df_m.to_excel(os.path.join(get_desktop_path(), raw_fname), index_label="Patient_ID")
    except: pass

    img_folder = os.path.join(get_desktop_path(), f"SubtypeCorr_{clean_filename(name1)}_{clean_filename(name2)}")
    if plot_enabled and not os.path.exists(img_folder): os.makedirs(img_folder)
    
    res = []
    groups = get_analysis_groups(df_m, subtype_col, subtypes_to_run)
            
    for grp_name, grp_df in groups:
        r, p, n = calculate_pearson(grp_df[name1], grp_df[name2])
        if r is not None:
            res.append({'Subtype':grp_name, 'N':n, 'R':r, 'P':p})
            if plot_enabled: 
                save_subtype_scatter(grp_df, name1, name2, img_folder, grp_name, r, p, n)
    
    if res:
        df_res = pd.DataFrame(res)
        if FDR_ENABLED:
            p_values = df_res['P'].dropna()
            if not p_values.empty:
                df_res['Adjusted_P_Value'] = multipletests(p_values, method='fdr_bh')[1]
        df_res.to_excel(os.path.join(get_desktop_path(), f"SubCorr_{clean_filename(name1)}_{clean_filename(name2)}.xlsx"), index=False)

def run_cox_analysis(target_name, target_data, df_clinical, subtype_col, time_col, status_col, plot_individual=True, save_excel=True, subtypes_to_run=None):
    merged_df = df_clinical.copy()
    merged_df['Target_Expression'] = pd.to_numeric(target_data, errors='coerce')
    if time_col not in merged_df.columns or status_col not in merged_df.columns:
        print("Error: Survival columns not found."); return []
    groups = get_analysis_groups(merged_df, subtype_col, subtypes_to_run)
    cox_results = []
    for grp_name, grp_df in groups:
        sub_df = grp_df[[time_col, status_col, 'Target_Expression']].dropna()
        if len(sub_df) < 2 or sub_df['Target_Expression'].std() == 0 or sub_df[status_col].nunique() < 2: continue
        
        # ---------------------------------------------------------
        # [Upgrade] Z-score scaling for continuous variables in Cox regression.
        # This ensures that the HR represents the risk increase per 1-SD change,
        # allowing for fair comparison between genes with different expression ranges in the forest plot.
        # ---------------------------------------------------------
        sub_df['Target_Expression'] = (sub_df['Target_Expression'] - sub_df['Target_Expression'].mean()) / sub_df['Target_Expression'].std()
        
        try:
            cph = CoxPHFitter(penalizer=0.05)
            cph.fit(sub_df, duration_col=time_col, event_col=status_col)
            summary = cph.summary.iloc[0]
            hr_val = cph.hazard_ratios_.iloc[0]
            if hr_val > 100 or hr_val < 0.01: continue
            cox_results.append({
                'Scope': grp_name, 'Gene/Signature': target_name,
                'Hazard_Ratio': hr_val,
                'CI_Lower_95': cph.confidence_intervals_.iloc[0, 0],
                'CI_Upper_95': cph.confidence_intervals_.iloc[0, 1],
                'P_Value': summary['p'],
                'N': len(sub_df), 'Events': sub_df[status_col].sum(),
                'Time_Col': time_col, 'Status_Col': status_col
            })
        except: pass
        
    if cox_results:
        if save_excel:
            res_df = pd.DataFrame(cox_results)
            fname = f"Cox_Hazard_{clean_filename(target_name)}_{clean_filename(time_col)}_{clean_filename(status_col)}.xlsx"
            res_df.to_excel(os.path.join(get_desktop_path(), fname), index=False)
        if plot_individual:
            draw_single_forest_plot(cox_results, target_name, time_col, status_col)
    return cox_results

def draw_single_forest_plot(cox_results, target_name, time_col='', status_col=''):
    if not cox_results: return
    setup_plot_style()
    df = pd.DataFrame(cox_results)
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['Hazard_Ratio', 'CI_Lower_95', 'CI_Upper_95'])
    if df.empty: return
    df = df[::-1].reset_index(drop=True)
    plt.figure(figsize=(9, 5 + len(df) * 0.8))
    x_err_lower = (df['Hazard_Ratio'] - df['CI_Lower_95']).abs()
    x_err_upper = (df['CI_Upper_95'] - df['Hazard_Ratio']).abs()
    plt.errorbar(df['Hazard_Ratio'], range(len(df)), xerr=[x_err_lower, x_err_upper], 
                 fmt='s', color='black', ecolor='gray', capsize=5, markersize=10)
    plt.axvline(x=1, color='red', linestyle='--', linewidth=1)
    plt.yticks(range(len(df)), df['Scope'], fontsize=18)
    for i, row in df.iterrows():
        p_val = row['P_Value']
        p_color = '#d62728' if p_val < 0.05 else 'black'
        fw = 'bold' if p_val < 0.05 else 'normal'
        p_txt = f"p {format_pval(p_val, use_e=False, with_stars=True)}"
        val_txt = f"{row['Hazard_Ratio']:.2f} [{row['CI_Lower_95']:.2f}-{row['CI_Upper_95']:.2f}]"
        plt.text(df['CI_Upper_95'].max() * 1.1, i, f"{val_txt}\n{p_txt}", verticalalignment='center', fontsize=15, color=p_color, fontweight=fw)
    title_suffix = f" ({time_col})" if time_col else ""
    plt.title(f"Cox Hazard Ratio: {target_name}{title_suffix}", fontsize=22, fontweight='bold')
    plt.xlabel("Hazard Ratio (log scale)", fontsize=20); plt.xscale('log')
    plt.tick_params(axis='x', which='major', labelsize=16)
    plt.xlim(max(0.1, df['CI_Lower_95'].min()*0.8), max(df['CI_Upper_95'].max()*1.5, 2.0))
    plt.tight_layout()
    file_suffix = f"_{clean_filename(time_col)}_{clean_filename(status_col)}" if time_col and status_col else ""
    plt.savefig(os.path.join(get_desktop_path(), f"Forest_Plot_{clean_filename(target_name)}{file_suffix}.{IMG_EXT}"), dpi=IMG_DPI)
    plt.close()

def draw_batch_forest_plot(batch_results, scope_name, orientation='1', time_col='', status_col=''):
    if not batch_results: return
    setup_plot_style()
    summary_data = []
    for res_list in batch_results:
        for item in res_list:
            if item['Scope'] == scope_name: summary_data.append(item)
    if not summary_data: return
    df = pd.DataFrame(summary_data)
    df = df.sort_values(by='Hazard_Ratio', ascending=True).reset_index(drop=True)
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['Hazard_Ratio', 'CI_Lower_95', 'CI_Upper_95'])
    if df.empty: return
    x_err_lower = (df['Hazard_Ratio'] - df['CI_Lower_95']).abs()
    x_err_upper = (df['CI_Upper_95'] - df['Hazard_Ratio']).abs()
    colors = ['#d62728' if hr > 1 else '#1f77b4' for hr in df['Hazard_Ratio']]
    
    if orientation == '2':
        plt.figure(figsize=(max(10, len(df) * 0.8), 6))
        for i in range(len(df)):
            plt.errorbar(i, df.loc[i, 'Hazard_Ratio'], yerr=[[x_err_lower[i]], [x_err_upper[i]]], 
                         fmt='o', color=colors[i], ecolor='gray', capsize=5, markersize=9)
        plt.axhline(y=1, color='black', linestyle='--', linewidth=1)
        plt.xticks(range(len(df)), df['Gene/Signature'], fontsize=16, fontweight='bold', rotation=45, ha='right')
        for i, row in df.iterrows():
            p_val = row['P_Value']
            p_color = '#d62728' if p_val < 0.05 else colors[i]
            fw = 'bold' if p_val < 0.05 else 'normal'
            p_txt = f"p {format_pval(p_val, use_e=False, with_stars=True)}"
            val_txt = f"{row['Hazard_Ratio']:.2f}"
            plt.text(i + 0.15, row['Hazard_Ratio'], f"{val_txt}\n({p_txt})", verticalalignment='center', fontsize=14, color=p_color, fontweight=fw)
        title_suffix = f" ({time_col})" if time_col else ""
        plt.title(f"Multi-Gene Survival Analysis ({scope_name}){title_suffix}", fontsize=22, fontweight='bold')
        plt.ylabel("Hazard Ratio (log scale)", fontsize=20); plt.yscale('log')
        plt.tick_params(axis='y', which='major', labelsize=16)
        plt.ylim(max(0.1, df['CI_Lower_95'].min() * 0.8), max(df['CI_Upper_95'].max() * 1.2, 2.0))
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        plt.tight_layout()
    else:
        plt.figure(figsize=(9, 5 + len(df) * 0.8))
        for i in range(len(df)):
            p_val = df.loc[i, 'P_Value']
            hr = df.loc[i, 'Hazard_Ratio']
            p_color = '#d62728' if p_val < 0.05 else colors[i]
            fw = 'bold' if p_val < 0.05 else 'normal'
            plt.errorbar(hr, i, xerr=[[x_err_lower[i]], [x_err_upper[i]]], 
                         fmt='o', color=colors[i], ecolor='gray', capsize=5, markersize=10)
            p_txt = f"p {format_pval(p_val, use_e=False, with_stars=True)}"
            val_txt = f"{hr:.2f}"
            plt.text(hr, i + 0.35, f"{val_txt} ({p_txt})", horizontalalignment='center', fontsize=14, color=p_color, fontweight=fw)
        plt.axvline(x=1, color='black', linestyle='--', linewidth=1)
        plt.yticks(range(len(df)), df['Gene/Signature'], fontsize=16, fontweight='bold')
        title_suffix = f" ({time_col})" if time_col else ""
        plt.title(f"Multi-Gene Survival Analysis ({scope_name}){title_suffix}", fontsize=22, fontweight='bold')
        plt.xlabel("Hazard Ratio (log scale)", fontsize=20); plt.xscale('log')
        plt.tick_params(axis='x', which='major', labelsize=16)
        plt.xlim(max(0.1, df['CI_Lower_95'].min() * 0.8), max(df['CI_Upper_95'].max() * 1.2, 2.0))
        plt.grid(axis='x', linestyle='--', alpha=0.5)
        plt.tight_layout()

    file_suffix = f"_{clean_filename(time_col)}_{clean_filename(status_col)}" if time_col and status_col else ""
    fname = f"Batch_Forest_{clean_filename(scope_name)}{file_suffix}.{IMG_EXT}"
    plt.savefig(os.path.join(get_desktop_path(), fname), dpi=IMG_DPI)
    plt.close()

def run_batch_correlation(list1, list2, df_gene1, df_path1, df_gene2, df_path2, df_clinical, subtype_col, suffix1="", suffix2="", subtypes_to_run=None, cmap='vlag'):
    print(f"\nStarting batch correlation analysis ({len(list1)} x {len(list2)})...")
    data_dict1 = {}
    data_dict2 = {}
    valid1, valid2 = [] , []
    for name in list1:
        d, t = get_data_by_name(name, df_gene1, df_path1)
        if d is not None: 
            display_name = name + suffix1 if suffix1 and 'Protein' not in name else name
            data_dict1[display_name] = d
            if display_name not in valid1: valid1.append(display_name)
    for name in list2:
        d, t = get_data_by_name(name, df_gene2, df_path2)
        if d is not None: 
            display_name = name + suffix2 if suffix2 and 'Protein' not in name else name
            data_dict2[display_name] = d
            if display_name not in valid2: valid2.append(display_name)
            
    if not valid1 or not valid2: return

    groups = get_analysis_indices(df_clinical, df_gene1.index, subtype_col, subtypes_to_run)

    fname = "Batch_Correlation_Matrix.xlsx"
    try:
        with pd.ExcelWriter(os.path.join(get_desktop_path(), fname), engine='openpyxl') as writer:
            for grp_name, grp_index in groups:
                print(f"-> Calculating: {grp_name}")
                results = []
                for n1 in valid1:
                    for n2 in valid2:
                                idx = grp_index.intersection(data_dict1[n1].index).intersection(data_dict2[n2].index)
                                if len(idx) < 3: continue
                                r, p, n = calculate_pearson(data_dict1[n1].loc[idx], data_dict2[n2].loc[idx])
                                if r is not None:
                                    results.append({'Name1': n1, 'Name2': n2, 'R': r, 'P': p, 'N': n})
                if results:
                    df_res = pd.DataFrame(results)
                    if FDR_ENABLED:
                        not_na_mask = df_res['P'].notna()
                        if not_na_mask.any():
                            df_res.loc[not_na_mask, 'Adjusted_P_Value'] = multipletests(df_res.loc[not_na_mask, 'P'], method='fdr_bh')[1]
                    df_res.attrs['cmap'] = cmap # Attach colormap info to the DataFrame
                    df_res.to_excel(writer, sheet_name=grp_name[:31], index=False)
                    draw_correlation_heatmap(df_res, grp_name)
    except: pass

def draw_correlation_heatmap(df_res, grp_name):
    if df_res.empty: return
    setup_plot_style()
    pivot_r = df_res.pivot_table(index='Name1', columns='Name2', values='R')
    pivot_p = df_res.pivot(index='Name1', columns='Name2', values='P')
    
    annot_matrix = []
    for r_row, p_row in zip(pivot_r.values, pivot_p.values):
        row_annot = []
        for r, p in zip(r_row, p_row):
            if pd.isna(r): row_annot.append("")
            else:
                star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                row_annot.append(f"{r:.2f}\n{star}")
        annot_matrix.append(row_annot)
        
    # [v82.2 Fix] Limit max figure size to prevent exceeding Matplotlib's pixel limit with too many genes (e.g., 44x44), which would fail to save.
    fig_width = min(60.0, max(8.0, len(pivot_r.columns) * 1.5))
    fig_height = min(60.0, max(6.0, len(pivot_r.index) * 1.5))
    plt.figure(figsize=(fig_width, fig_height))
    heatmap_ax = sns.heatmap(pivot_r, annot=annot_matrix, fmt="", cmap=df_res.attrs.get('cmap', 'vlag'), vmin=-1, vmax=1,
                linewidths=0.5, linecolor='white', annot_kws={"size": 14})
    cbar = heatmap_ax.collections[0].colorbar
    if cbar:
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label('Correlation (R / Rho)', size=18)
    plt.title(f"Correlation Matrix ({grp_name})", fontsize=22, fontweight='bold', pad=20)
    plt.xticks(rotation=45, ha='right', fontsize=16, fontweight='bold')
    plt.yticks(rotation=0, fontsize=16, fontweight='bold')
    plt.xlabel(""); plt.ylabel("")
    plt.tight_layout()
    fname = f"Heatmap_Matrix_{clean_filename(grp_name)}.{IMG_EXT}"
    plt.savefig(os.path.join(get_desktop_path(), fname), dpi=IMG_DPI)
    plt.close()

def run_patient_heatmap(target_name, target_data, gene_list, df_gene, df_pathway, df_clinical, subtype_col, orientation, sort_order='1', subtypes_to_run=None, suffix="", hl_set=None, cmap='vlag'):
    print(f"\n>>> Preparing to plot patient co-expression heatmap for {target_name}...")
    data_dict = {}
    for g in gene_list:
        d, t = get_data_by_name(g, df_gene, df_pathway)
        if d is not None: 
            display_g = g + suffix if suffix else g
            data_dict[display_g] = d
    if not data_dict: return
        
    groups = get_analysis_indices(df_clinical, target_data.index, subtype_col, subtypes_to_run)
            
    for grp_name, grp_index in groups:
        valid_idx = grp_index.intersection(target_data.dropna().index)
        if len(valid_idx) < 2: continue
        sub_target = target_data.loc[valid_idx]
        
        corr_results = []
        for g, d in data_dict.items():
            g_data = d.loc[valid_idx]
            r, p, n = calculate_pearson(sub_target, g_data)
            if r is not None: corr_results.append({'Gene': g, 'R': r, 'P': p, 'Data': g_data})
                
        if not corr_results: continue
        is_reverse = True if sort_order == '1' else False
        corr_results.sort(key=lambda x: x['R'], reverse=is_reverse)

        # [User Request] Output correlation summary to Excel
        summary_df = pd.DataFrame([{'Gene': r['Gene'], 'R': r['R'], 'P': r['P']} for r in corr_results])
        if not summary_df.empty:
            if FDR_ENABLED:
                p_values = summary_df['P'].dropna()
                if not p_values.empty:
                    # Calculate FDR using multipletests
                    summary_df['FDR'] = multipletests(p_values, method='fdr_bh')[1]
            
            summary_fname = f"Summary_PatHeatmap_{clean_filename(target_name)}_{clean_filename(grp_name)}.xlsx"
            try:
                summary_df.to_excel(os.path.join(get_desktop_path(), summary_fname), index=False)
            except Exception as e:
                print(f"  -> [Warning] Could not save summary Excel file: {e}")
        
        boundary_pos = None
        for i in range(len(corr_results) - 1):
            r_curr = corr_results[i]['R']
            r_next = corr_results[i+1]['R']
            if (r_curr >= 0 > r_next) or (r_curr < 0 <= r_next):
                boundary_pos = i + 2 
                break

        sorted_patients = sub_target.sort_values(ascending=True).index
        plot_data = []
        row_labels = [f"{target_name} (Target)"]
        row_colors = []
        row_weights = []
        
        t_clean = target_name.replace(" (Protein)", "").strip()
        if hl_set and normalize_name(t_clean) in hl_set:
            row_colors.append('#D62728')
            row_weights.append('bold')
        else:
            row_colors.append('black')
            row_weights.append('normal')
        t_z = (sub_target - sub_target.mean()) / (sub_target.std() + 1e-9)
        plot_data.append(t_z.loc[sorted_patients].values)
        
        for item in corr_results:
            g = item['Gene']; r = item['R']; p = item['P']; g_data = item['Data']
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            row_labels.append(f"{g} (R={r:.2f} {star})")
            z = (g_data - g_data.mean()) / (g_data.std() + 1e-9)
            plot_data.append(z.loc[sorted_patients].values)
            
            g_clean = g.replace(suffix, '').strip() if suffix else g
            if hl_set and normalize_name(g_clean) in hl_set:
                row_colors.append('#D62728')
                row_weights.append('bold')
            else:
                row_colors.append('black')
                row_weights.append('normal')
            
        matrix = pd.DataFrame(plot_data, index=row_labels, columns=sorted_patients).fillna(0)
        
        try:
            raw_fname = f"RawData_PatHeatmap_{clean_filename(target_name)}_{clean_filename(grp_name)}.xlsx"
            matrix.to_excel(os.path.join(get_desktop_path(), raw_fname))
        except: pass
        
        setup_plot_style()
        if orientation == '2': 
            matrix_plot = matrix.T
            fig_w = max(6, len(row_labels) * 0.5)
            fig_h = max(8, min(15, len(sorted_patients) * 0.05)) 
            plt.figure(figsize=(fig_w, fig_h))
            
            ax = sns.heatmap(matrix_plot, cmap=cmap, vmin=-3, vmax=3,
                        xticklabels=True, yticklabels=False, 
                        cbar_kws={'label': 'Expression (Z)'})
            ax.invert_yaxis()
            ax.axvline(1, color='black', linewidth=3)
            
            if boundary_pos is not None:
                ax.axvline(boundary_pos, color='#333333', linewidth=2.5, linestyle='--')
                grp1_label = "Positive Corr (+)" if corr_results[0]['R'] >= 0 else "Negative Corr (-)"
                grp2_label = "Negative Corr (-)" if corr_results[0]['R'] >= 0 else "Positive Corr (+)"
                grp1_color = "#d62728" if corr_results[0]['R'] >= 0 else "#1f77b4"
                grp2_color = "#1f77b4" if corr_results[0]['R'] >= 0 else "#d62728"
                
                y_bracket = len(sorted_patients) * 1.02  
                
                ax.annotate('', xy=(1.05, y_bracket), xytext=(boundary_pos - 0.05, y_bracket),
                            xycoords='data', textcoords='data', annotation_clip=False,
                            arrowprops=dict(arrowstyle='-', color=grp1_color, lw=2.5))
                ax.text((1 + boundary_pos)/2, y_bracket + len(sorted_patients)*0.015, grp1_label, 
                        color=grp1_color, ha='center', va='bottom', fontweight='bold', fontsize=15)
                
                ax.annotate('', xy=(boundary_pos + 0.05, y_bracket), xytext=(len(row_labels) - 0.05, y_bracket),
                            xycoords='data', textcoords='data', annotation_clip=False,
                            arrowprops=dict(arrowstyle='-', color=grp2_color, lw=2.5))
                ax.text((boundary_pos + len(row_labels))/2, y_bracket + len(sorted_patients)*0.015, grp2_label,
                        color=grp2_color, ha='center', va='bottom', fontweight='bold', fontsize=15)

            plt.title(f"Patient Co-expression: {target_name}\n({grp_name}, N={len(sorted_patients)})", fontsize=22, fontweight='bold', pad=90)
            plt.xticks(rotation=45, ha='right', fontsize=16)
            for tick, color, weight in zip(ax.get_xticklabels(), row_colors, row_weights):
                tick.set_color(color)
                tick.set_fontweight(weight)
            plt.ylabel("Patients (Sorted by Target ->)", fontsize=18)
            
        else: 
            fig_w = max(10, min(18, len(sorted_patients) * 0.05))
            fig_h = max(6, len(row_labels) * 0.5)
            plt.figure(figsize=(fig_w, fig_h))
            
            ax = sns.heatmap(matrix, cmap=cmap, vmin=-3, vmax=3,
                        xticklabels=False, yticklabels=True, 
                        cbar_kws={'label': 'Expression (Z)', 'orientation': 'horizontal', 'pad': 0.12})
            ax.axhline(1, color='black', linewidth=3)
            
            if boundary_pos is not None:
                ax.axhline(boundary_pos, color='#333333', linewidth=2.5, linestyle='--')
                grp1_label = "Positive Corr (+)" if corr_results[0]['R'] >= 0 else "Negative Corr (-)"
                grp2_label = "Negative Corr (-)" if corr_results[0]['R'] >= 0 else "Positive Corr (+)"
                grp1_color = "#d62728" if corr_results[0]['R'] >= 0 else "#1f77b4"
                grp2_color = "#1f77b4" if corr_results[0]['R'] >= 0 else "#d62728"
                
                x_bracket = len(sorted_patients) * 1.05  
                
                ax.annotate('', xy=(x_bracket, 1.05), xytext=(x_bracket, boundary_pos - 0.05),
                            xycoords='data', textcoords='data', annotation_clip=False,
                            arrowprops=dict(arrowstyle='-', color=grp1_color, lw=2.5))
                ax.text(x_bracket + len(sorted_patients)*0.02, (1 + boundary_pos)/2, grp1_label,
                        color=grp1_color, ha='left', va='center', fontweight='bold', fontsize=15, rotation=-90)
                
                ax.annotate('', xy=(x_bracket, boundary_pos + 0.05), xytext=(x_bracket, len(row_labels) - 0.05),
                            xycoords='data', textcoords='data', annotation_clip=False,
                            arrowprops=dict(arrowstyle='-', color=grp2_color, lw=2.5))
                ax.text(x_bracket + len(sorted_patients)*0.02, (boundary_pos + len(row_labels))/2, grp2_label,
                        color=grp2_color, ha='left', va='center', fontweight='bold', fontsize=15, rotation=-90)

            plt.title(f"Patient Co-expression: {target_name}\n({grp_name}, N={len(sorted_patients)})", fontsize=22, fontweight='bold', pad=50)
            plt.yticks(rotation=0, fontsize=16)
            for tick, color, weight in zip(ax.get_yticklabels(), row_colors, row_weights):
                tick.set_color(color)
                tick.set_fontweight(weight)
            plt.xlabel("Patients (Sorted by Target ->)", fontsize=18)
            
        plt.tight_layout() 
        img_fname = f"PatHeatmap_{clean_filename(target_name)}_{clean_filename(grp_name)}.{IMG_EXT}"
        plt.savefig(os.path.join(get_desktop_path(), img_fname), dpi=IMG_DPI, bbox_inches='tight') 
        plt.close()

# [v81.1 修正] Matplotlib 強制複數參數：linewidths & edgecolors
def draw_ext_bar_lollipop(ax, df, title, highlight_set, color_high, color_norm, chart_style, val_col, val_label, x_col='mlog10p', x_label='-log10(P-value)', label_pad=5, size_col=None, size_range=(200, 800)):
    df = df.copy()
    y_pos = np.arange(len(df))
    colors = [color_high if normalize_name(tf) in highlight_set else color_norm for tf in df['Item_Name']]
    df['is_high'] = df['Item_Name'].apply(lambda tf: normalize_name(tf) in highlight_set)
    
    x_max = df[x_col].max()
    x_min = df[x_col].min()
    data_range = x_max - x_min if x_max != x_min else 1
    
    if chart_style == "bar": 
        ax.barh(y_pos, df[x_col], color=colors, edgecolor='black', linewidth=0.8, height=0.7, zorder=2)
    else:
        ax.hlines(y=y_pos, xmin=0, xmax=df[x_col], color=colors, alpha=0.6, linewidth=2, zorder=1)
        if size_col and size_col in df.columns:
            # Fix: change edgecolor to edgecolors, linewidth to linewidths
            sns.scatterplot(data=df, x=x_col, y=y_pos, hue='is_high', size=size_col, 
                            sizes=size_range, palette={True: color_high, False: color_norm},
                            edgecolors='black', linewidths=1, ax=ax, legend='brief', zorder=3)
            
            handles, labels = ax.get_legend_handles_labels()
            if size_col in labels:
                size_idx = labels.index(size_col)
                size_handles = handles[size_idx+1:]
                raw_size_labels = labels[size_idx+1:]
                
                clean_size_labels = []
                for l in raw_size_labels:
                    try:
                        lf = float(l)
                        if abs(lf - round(lf)) < 1e-6:
                            clean_size_labels.append(str(int(round(lf))))
                        else:
                            clean_size_labels.append(f"{lf:.2f}")
                    except:
                        clean_size_labels.append(l)
                        
                if ax.get_legend() is not None:
                    ax.get_legend().remove()
                
                legend = ax.legend(size_handles, clean_size_labels, title=size_col.replace('_', ' '), loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False, labelspacing=1.2)
                for handle in legend.legend_handles:
                    if hasattr(handle, 'set_facecolors'):
                        handle.set_facecolors('#E0E0E0')
                        handle.set_edgecolors('black')
            else:
                if ax.get_legend() is not None:
                    ax.get_legend().remove()
        else:
            # Same fix here: use edgecolors and linewidths
            ax.scatter(df[x_col], y_pos, color=colors, s=250, edgecolors='black', linewidths=1, zorder=3)
        
    if x_col == 'mlog10p':
        sig_line = -np.log10(0.05)
        ax.axvline(x=sig_line, color='#D62728', linestyle='--', alpha=0.8, linewidth=1.5, zorder=5)
        ax.text(sig_line, ax.get_ylim()[1], "p=0.05", color='#D62728', fontsize=13, fontweight='bold', ha='center', va='bottom', zorder=6)
    else:
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.3, linewidth=1.5, zorder=1)

    ax.set_xlabel(x_label, fontsize=16, fontweight='bold')
    ax.set_yticks(y_pos)
    display_names = [str(name).replace('_', ' ') for name in df['Item_Name']]
    ax.set_yticklabels(display_names, fontsize=18, fontweight='bold')
    for tick in ax.get_yticklabels():
        if normalize_name(tick.get_text()) in highlight_set: tick.set_color(color_high)
    ax.tick_params(axis='y', pad=label_pad) 
    ax.grid(axis='x', linestyle=':', alpha=0.3, zorder=0)
    
    if chart_style != "bar":
        ax.set_xlim(min(0, x_min - data_range * 0.05), max(0, x_max + data_range * 0.6))
    else:
        ax.set_xlim(min(0, x_min - data_range * 0.05), max(0, x_max + data_range * 0.3))

    for i, (_, row) in enumerate(df.iterrows()):
        is_high = row['is_high']
        val_str = ""
        if val_col and val_col in row and pd.notna(row[val_col]):
            val = row[val_col]
            try: 
                val_f = float(val)
                val_col_lower = val_col.lower().strip()
                
                is_pvalue = val_col_lower in ['p', 'pvalue', 'p-value', 'fdr', 'padj', 'qvalue', 'q-value']
                
                if is_pvalue:
                    if val_f < 0.0001:
                        val_formatted = "< 0.0001"
                    else:
                        val_formatted = f"{val_f:.3e}"
                else:
                    if abs(val_f - round(val_f)) < 1e-6:
                        val_formatted = str(int(round(val_f)))
                    else:
                        val_formatted = f"{val_f:.2f}"
                val_str = f"{val_label}: {val_formatted}" if val_label else val_formatted
            except: 
                val_str = f"{val_label}: {val}" if val_label else f"{val}"
            
        if chart_style != "bar" and size_col and size_col in df.columns:
            s_val = row[size_col]
            s_max_val = df[size_col].max()
            if pd.notna(s_val) and s_max_val > 0:
                radius_frac = np.sqrt(float(s_val)) / np.sqrt(float(s_max_val))
                row_offset = data_range * 0.04 + (data_range * 0.12) * radius_frac 
            else:
                row_offset = data_range * 0.16
        elif chart_style != "bar":
            row_offset = data_range * 0.12
        else:
            row_offset = data_range * 0.02
                
        text_x = row[x_col] + row_offset if row[x_col] >= 0 else row[x_col] - row_offset
        ha = 'left' if row[x_col] >= 0 else 'right'
        if val_str:
            ax.text(text_x, i, val_str, va='center', ha=ha, fontsize=15, fontweight='bold' if is_high else 'normal', color=color_high if is_high else '#333333', zorder=7)
            
    if title:
        ax.set_title(title, fontsize=24, fontweight='bold', pad=35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

def draw_external_waterfall(ax, df_plot, val_col, title, hl_set, c_pos, c_neg, is_horizontal, user_min=None, user_max=None):
    """
    Core drawing function for the external waterfall plot.
    """
    pos = np.arange(len(df_plot))
    colors = []
    for v, name in zip(df_plot[val_col], df_plot['Item_Name']):
        if hl_set:
            if normalize_name(name) in hl_set:
                colors.append(c_pos if v >= 0 else c_neg)
            else:
                colors.append('#E0E0E0')
        else:
            colors.append(c_pos if v >= 0 else c_neg)
    
    bar_lw = 0.8 if len(df_plot) <= 60 else (0.2 if len(df_plot) <= 150 else 0)
    e_color = 'black' if bar_lw > 0 else 'none'

    if is_horizontal:
        bars = ax.barh(pos, df_plot[val_col], color=colors, edgecolor=e_color, linewidth=bar_lw, height=1.0)
    else:
        bars = ax.bar(pos, df_plot[val_col], color=colors, edgecolor=e_color, linewidth=bar_lw, width=1.0)
    
    val_max = df_plot[val_col].max() if not df_plot.empty else 0
    val_min = df_plot[val_col].min() if not df_plot.empty else 0
    val_range = val_max - val_min if val_max != val_min else 1.0

    # Determine axis limits
    if user_min is None:
        expand_factor = 1.6 if hl_set else 0.2
        user_min = float(val_min - val_range * expand_factor)
    if user_max is None:
        expand_factor = 1.6 if hl_set else 0.2
        user_max = float(val_max + val_range * expand_factor)

    if hl_set:
        offset = val_range * 0.05

        if is_horizontal: ax.set_yticks([])
        else: ax.set_xticks([])

        texts = []
        target_xs = []
        target_ys = []
        n_bars = len(df_plot)
        import textwrap
        n_hl = sum(1 for name in df_plot['Item_Name'] if normalize_name(name) in hl_set)
        fs = 13 if n_hl <= 15 else (11 if n_hl <= 30 else (9 if n_hl <= 60 else 7))
        
        stagger_levels = min(20, max(4, n_hl // 3))
        avail_space_pos = user_max - val_max if user_max > val_max else val_range * 0.2
        avail_space_neg = val_min - user_min if val_min > user_min else val_range * 0.2

        for i, (v, name) in enumerate(zip(df_plot[val_col], df_plot['Item_Name'])):
            if normalize_name(name) in hl_set:
                clean_n = str(name).replace('_', ' ')
                clean_n = '\n'.join(textwrap.wrap(clean_n, width=14))
                t_c = c_pos if v >= 0 else c_neg
                
                is_pos_placement = len(texts) % 2 == 0
                if is_pos_placement:
                    stagger_step = (avail_space_pos * 0.85) / stagger_levels if n_hl > 15 else (val_range * 0.12)
                else:
                    stagger_step = (avail_space_neg * 0.85) / stagger_levels if n_hl > 15 else (val_range * 0.12)
                    
                stagger = ((len(texts) // 2) % stagger_levels) * stagger_step + np.random.uniform(-0.01, 0.01) * val_range
                
                if is_horizontal:
                    if i < n_bars * 0.25: va = 'bottom'
                    elif i > n_bars * 0.75: va = 'top'
                    else: va = 'center'
                    if is_pos_placement:
                        x_pos_text = val_max + offset + stagger
                    else:
                        x_pos_text = val_min - offset - stagger
                    target_x = v if abs(x_pos_text - v) < abs(x_pos_text - 0) else 0
                    t = ax.text(x_pos_text, i, clean_n, ha='left' if is_pos_placement else 'right', va=va, rotation=0, fontsize=fs, fontweight='900', color=t_c, zorder=10, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
                    texts.append(t)
                    target_xs.append(target_x)
                    target_ys.append(i)
                else:
                    if i < n_bars * 0.25: ha = 'left'
                    elif i > n_bars * 0.75: ha = 'right'
                    else: ha = 'center'
                    if is_pos_placement:
                        y_pos_text = val_max + offset + stagger
                    else:
                        y_pos_text = val_min - offset - stagger
                    target_y = v if abs(y_pos_text - v) < abs(y_pos_text - 0) else 0
                    t = ax.text(i, y_pos_text, clean_n, ha=ha, va='bottom' if is_pos_placement else 'top', rotation=0, fontsize=fs, fontweight='900', color=t_c, zorder=10, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
                    texts.append(t)
                    target_xs.append(i)
                    target_ys.append(target_y)

        if texts:
            try: ax.figure.canvas.draw()
            except Exception: pass
            
            if HAS_ADJUST_TEXT:
                import contextlib
                with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                    try:
                        # 移除 x=target_xs, y=target_ys，避免 adjust_text 將標籤強拉回中心點
                        adjust_text(texts, force_text=(0.3, 0.3), expand_text=(1.1, 1.1), max_iter=5000)
                    except Exception:
                        pass
            else:
                print("❌ Warning: adjustText not installed! Labels may overlap. Please run: pip install adjustText")

            for t, tx, ty in zip(texts, target_xs, target_ys):
                ax.annotate("", xy=(tx, ty), xytext=t.get_position(),
                            xycoords='data', textcoords='data',
                            arrowprops=dict(arrowstyle='-', color=t.get_color(), lw=1.2, alpha=0.75, patchA=t, shrinkB=0),
                            zorder=5)
    else:
        display_names = [str(n).replace('_',' ') for n in df_plot['Item_Name']]
        if len(df_plot) > 60:
            if is_horizontal: ax.set_yticks([])
            else: ax.set_xticks([])
        else:
            if is_horizontal:
                ax.set_yticks(pos)
                ax.set_yticklabels(display_names, ha='right', fontsize=13)
            else:
                ax.set_xticks(pos)
                ax.set_xticklabels(display_names, rotation=90, ha='center', fontsize=13)
    
    if is_horizontal:
        ax.set_ylim(-0.5, len(df_plot) - 0.5)
        ax.set_xlabel(val_col, fontsize=18, fontweight='bold')
        ax.spines['left'].set_visible(False)
        ax.set_xlim(user_min, user_max)
    else:
        ax.set_xlim(-0.5, len(df_plot) - 0.5)
        ax.set_ylabel(val_col, fontsize=18, fontweight='bold')
        ax.set_ylim(user_min, user_max)

    if title:
        ax.set_title(title, fontsize=22, fontweight='bold', pad=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def draw_ext_volcano(ax, df, title, target_genes, y_thresh, y_dir, x_thresh_r, x_thresh_l, x_name, y_name, super_targets=None):
    if str(y_dir) == '1':
        sig_up = (df['Plot_Y'] >= y_thresh) & (df['Plot_X'] >= x_thresh_r)
        sig_down = (df['Plot_Y'] >= y_thresh) & (df['Plot_X'] <= x_thresh_l)
    else:
        sig_up = (df['Plot_Y'] <= y_thresh) & (df['Plot_X'] >= x_thresh_r)
        sig_down = (df['Plot_Y'] <= y_thresh) & (df['Plot_X'] <= x_thresh_l)
        
    not_sig = ~(sig_up | sig_down)
    
    ax.scatter(df.loc[not_sig, 'Plot_X'], df.loc[not_sig, 'Plot_Y'], color='grey', alpha=0.4, s=15, label='Not Sig')
    ax.scatter(df.loc[sig_down, 'Plot_X'], df.loc[sig_down, 'Plot_Y'], color='#1F77B4', alpha=0.7, s=30, label='Down')
    ax.scatter(df.loc[sig_up, 'Plot_X'], df.loc[sig_up, 'Plot_Y'], color='#D62728', alpha=0.7, s=30, label='Up')
    
    ax.axhline(y=y_thresh, color='black', linestyle='--', linewidth=1, alpha=0.6, zorder=1)
    ax.axvline(x=x_thresh_r, color='black', linestyle='--', linewidth=1, alpha=0.6, zorder=1)
    ax.axvline(x=x_thresh_l, color='black', linestyle='--', linewidth=1, alpha=0.6, zorder=1)
    
    df['Norm_Name'] = df['Item_Name'].apply(normalize_name)
    super_norm = set([normalize_name(t) for t in (super_targets or [])])
    
    # Calculate data range to use as a basis for small perturbations
    x_range = df['Plot_X'].max() - df['Plot_X'].min()
    y_range = df['Plot_Y'].max() - df['Plot_Y'].min()
    x_range = x_range if x_range != 0 else 1
    y_range = y_range if y_range != 0 else 1

    
    # Auto-adjust font size and repulsion range based on the number of labels
    n_labels = sum(df['Norm_Name'].isin(super_norm | set([normalize_name(t) for t in target_genes])))
    base_fontsize = 16
    if n_labels >= 50: base_fontsize = 10
    elif n_labels >= 30: base_fontsize = 12
    elif n_labels >= 15: base_fontsize = 14

    texts_super = []
    texts_norm = []

    for tf in target_genes:
        target_norm = normalize_name(tf)
        gene_data = df[df['Norm_Name'] == target_norm]
        for _, row in gene_data.iterrows():
            display_name = str(row['Item_Name']).replace('_', ' ')
            
            px = row['Plot_X']
            py = row['Plot_Y']
            
            if target_norm in super_norm:
                # Strong highlight style (gold background, red border, dark red text, plus a halo)
                t = ax.text(px, py, display_name,
                            fontsize=base_fontsize+2, fontweight='900', color='#8B0000', 
                            bbox=dict(boxstyle="round,pad=0.3", fc="#FFFFCC", ec="#D62728", lw=1.5, alpha=0.95), 
                            zorder=15)
                ax.scatter(px, py, color='#FFFF00', edgecolors='#D62728', linewidths=2.5, s=150, zorder=14)
                texts_super.append(t)
            else:
                # Normal highlight style (white background, gray border, black text)
                ax.scatter(px, py, color='#39FF14', edgecolors='black', linewidths=1.5, s=80, zorder=13)
                t = ax.text(px, py, display_name,
                            fontsize=base_fontsize, fontweight='bold', color='black', 
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.8), 
                            zorder=10)
                texts_norm.append(t)
            
    if texts_super or texts_norm:
        c_xmin, c_xmax = ax.get_xlim()
        c_ymin, c_ymax = ax.get_ylim()
        x_expand = 0.25 if n_labels > 20 else 0.1
        y_expand_up = 0.25 if n_labels > 20 else 0.1
        y_expand_dn = 0.05
        ax.set_xlim(c_xmin - x_range * x_expand, c_xmax + x_range * x_expand)
        ax.set_ylim(c_ymin - y_range * y_expand_dn, c_ymax + y_range * y_expand_up)

    if HAS_ADJUST_TEXT and (texts_super or texts_norm): 
        # Force Matplotlib to render first to ensure the coordinate transformation matrix is initialized
        try:
            ax.figure.canvas.draw()
        except Exception:
            pass
            
        import contextlib
        # Create a silent block to suppress stubborn warnings from adjustText
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
            
            # Phase 1: First, lay out the super-highlighted genes (let them take the best positions closest to their points)
            if texts_super:
                adjust_text(texts_super, ax=ax, 
                            arrowprops=dict(arrowstyle='-', color='#D62728', lw=1.5, alpha=0.9),
                            force_text=(0.5, 1.0),
                            force_points=(1.0, 1.5),
                            expand_text=(1.2, 1.2),
                            expand_points=(1.2, 1.2),
                            max_iter=3000)
            
            # Phase 2: Lay out normal genes, setting super-highlighted ones as immovable obstacles (add_objects)
            if texts_norm:
                f_txt = (2.5, 4.0) if n_labels > 20 else (0.8, 1.2)
                f_pts = (3.0, 5.0) if n_labels > 20 else (1.0, 1.5)
                e_txt = (1.5, 1.8) if n_labels > 20 else (1.2, 1.3)
                adjust_text(texts_norm, ax=ax, 
                            add_objects=texts_super,  # The core magic: prevent normal labels from overlapping super labels
                            arrowprops=dict(arrowstyle='-', color='#777777', lw=1.0, alpha=0.8),
                            force_text=f_txt,
                            force_points=f_pts,
                            expand_text=e_txt,
                            expand_points=(1.5, 1.5),
                            max_iter=8000)
    
    ax.set_xlabel(x_name, fontsize=20, fontweight='bold')
    ax.set_ylabel(y_name, fontsize=20, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=24, fontweight='bold', pad=55) 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False, fontsize=16)

def run_correlation_screening(target_name, target_data, df_gene, df_pathway, p_val_thresh=0.05, save_enabled=True, gene_to_drugs=None, pathway_to_genes=None):
    print(f"\nRunning genome-wide correlation screening for {target_name} (P < {p_val_thresh}) (All Patients)...")
    def _scan_correlations(target_data, df_pool, desc):
        res = []
        target_data = pd.to_numeric(target_data, errors='coerce')
        for col_name, col_data in tqdm(df_pool.items(), total=df_pool.shape[1], desc=desc, ncols=80, leave=False):
            r, p, n = calculate_pearson(target_data, col_data)
            if r is not None: res.append({'Name': col_name, 'R': r, 'P': p, 'N': n})
        return pd.DataFrame(res)

    def get_top_tables(df_res, type_name):
        if df_res.empty: return pd.DataFrame(), pd.DataFrame()
        df_sig = df_res[df_res['P'] < p_val_thresh].copy()        
        # [v82.2 Fix] Add a hard limit to prevent memory explosion if p-value threshold is too loose.
        pos = df_sig[df_sig['R'] > 0].sort_values(by='R', ascending=False)
        neg = df_sig[df_sig['R'] < 0].sort_values(by='R', ascending=True)
        pos = pos.head(5000)
        neg = neg.head(5000)

        pos_final = pos.rename(columns={'Name': f'Pos_{type_name}_Name', 'R': f'Pos_{type_name}_R', 'P': f'Pos_{type_name}_P', 'N': f'Pos_{type_name}_N', 'FDR': f'Pos_{type_name}_FDR'})
        neg_final = neg.rename(columns={'Name': f'Neg_{type_name}_Name', 'R': f'Neg_{type_name}_R', 'P': f'Neg_{type_name}_P', 'N': f'Neg_{type_name}_N', 'FDR': f'Neg_{type_name}_FDR'})
        return pos_final.reset_index(drop=True), neg_final.reset_index(drop=True)

    df_g = _scan_correlations(target_data, df_gene, "Scanning Genes")
    df_p = _scan_correlations(target_data, df_pathway, "Scanning Paths")

    # [FDR] Combine p-values from genes and pathways for FDR correction
    all_results = pd.concat([df_g, df_p], ignore_index=True)
    if FDR_ENABLED:
        if not all_results.empty:
            all_results['FDR'] = multipletests(all_results['P'].dropna(), method='fdr_bh')[1] if not all_results['P'].dropna().empty else np.nan
            df_g = all_results.iloc[:len(df_g)].copy()
            df_p = all_results.iloc[len(df_g):].copy()

    pg_pos, pg_neg = get_top_tables(df_g, "Gene")
    pp_pos, pp_neg = get_top_tables(df_p, "Path")
    
    pg_pos = annotate_top_tables(pg_pos, "Gene", "Pos_Gene", gene_to_drugs, pathway_to_genes) # FDR is already in the df
    pg_neg = annotate_top_tables(pg_neg, "Gene", "Neg_Gene", gene_to_drugs, pathway_to_genes)
    pp_pos = annotate_top_tables(pp_pos, "Path", "Pos_Path", gene_to_drugs, pathway_to_genes)
    pp_neg = annotate_top_tables(pp_neg, "Path", "Neg_Path", gene_to_drugs, pathway_to_genes)
    if not pg_pos.empty: print(pg_pos.iloc[:, 0:2].head(5).to_string(index=False))
    if save_enabled:
        fname = f"Screening_{clean_filename(target_name)}.xlsx"
        try:
            with pd.ExcelWriter(os.path.join(get_desktop_path(), fname), engine='openpyxl') as writer:
                final = pd.concat([pg_pos, pg_neg, pp_pos, pp_neg], axis=1)
                final.to_excel(writer, sheet_name='Sig_Correlation', index=False)
        except: pass

def run_subtype_screening(target_name, target_data, df_gene, df_pathway, df_clinical, subtype_col, p_val_thresh=0.05, subtypes_to_run=None, gene_to_drugs=None, pathway_to_genes=None):    
    # [v82.1 Fix] If no stratification (global or local) is selected by the user, automatically switch to the global screening mode (run_correlation_screening).
    # This prevents an ExcelWriter error "At least one sheet must be visible" when there are no sheets to write.
    is_no_subtype_run = False
    if subtype_col is None:
        if subtypes_to_run is None or (len(subtypes_to_run) == 1 and subtypes_to_run[0] is None):
            is_no_subtype_run = True

    if is_no_subtype_run:
        run_correlation_screening(target_name, target_data, df_gene, df_pathway, p_val_thresh=p_val_thresh, save_enabled=True, gene_to_drugs=gene_to_drugs, pathway_to_genes=pathway_to_genes) # Fallback to global screening
        return
        
    print(f"\n【Stratified Screening】Analyzing {target_name} (P < {p_val_thresh})...")

    base_index = target_data.index
    if df_gene is not None: base_index = base_index.intersection(df_gene.index)
    if df_pathway is not None: base_index = base_index.intersection(df_pathway.index)

    groups = get_analysis_indices(df_clinical, base_index, subtype_col, subtypes_to_run)
    fname = f"Subtype_Screening_{clean_filename(target_name)}.xlsx"
    sheet_written = False  # Add a flag to track if any sheet has been written
    try:
        with pd.ExcelWriter(os.path.join(get_desktop_path(), fname), engine='openpyxl') as writer:
            for grp_name, grp_index in groups:                
                try:                    
                    print(f"--> Calculating: {grp_name}")                    
                    # [v82.2 Fix] Ensure all data is aligned based on the samples in the current group
                    valid_samples = grp_index.intersection(target_data.index).intersection(df_gene.index).intersection(df_pathway.index)
                    if len(valid_samples) < 3: continue

                    sub_target = target_data.loc[valid_samples]
                    sub_gene = df_gene.loc[valid_samples] if df_gene is not None else pd.DataFrame()
                    sub_path = df_pathway.loc[valid_samples] if df_pathway is not None else pd.DataFrame()
                    
                    res_g = []
                    target_num = pd.to_numeric(sub_target, errors='coerce')
                    for col_name, col_data in tqdm(sub_gene.items(), total=sub_gene.shape[1], desc="Genes", ncols=80, leave=False):
                        r, p, n = calculate_pearson(target_num, col_data)
                        if r is not None: res_g.append({'Name': col_name, 'R': r, 'P': p, 'N': n})
                    
                    res_p = []
                    for col_name, col_data in tqdm(sub_path.items(), total=sub_path.shape[1], desc="Paths", ncols=80, leave=False):
                        r, p, n = calculate_pearson(target_num, col_data)
                        if r is not None: res_p.append({'Name': col_name, 'R': r, 'P': p, 'N': n})
                    
                    df_g = pd.DataFrame(res_g)
                    df_p = pd.DataFrame(res_p)
                    if FDR_ENABLED:
                        all_results = pd.concat([df_g, df_p], ignore_index=True)
                        if not all_results.empty:
                            # [v82.3 Fix] Handle NaN P-values to avoid length mismatch after FDR correction
                            not_na_mask = all_results['P'].notna()
                            if not_na_mask.any():
                                p_values_no_na = all_results.loc[not_na_mask, 'P']
                                fdr_values = multipletests(p_values_no_na, method='fdr_bh')[1]
                                all_results['Adjusted_P_Value'] = np.nan
                                all_results.loc[not_na_mask, 'Adjusted_P_Value'] = fdr_values
                        df_g = all_results.iloc[:len(df_g)].copy()
                        df_p = all_results.iloc[len(df_g):].copy()

                    def get_top_tables(df_res, type_name): # FDR is now Adjusted_P_Value
                        if df_res.empty: return pd.DataFrame(), pd.DataFrame()
                        df_sig = df_res[df_res['P'] < p_val_thresh].copy()
                        # [v82.2 Fix] Add a hard limit to prevent memory explosion if p-value threshold is too loose.
                        pos = df_sig[df_sig['R'] > 0].sort_values(by='R', ascending=False).head(5000)
                        neg = df_sig[df_sig['R'] < 0].sort_values(by='R', ascending=True).head(5000)
                        
                        pos_final = pos.rename(columns={'Name': f'Pos_{type_name}_Name', 'R': f'Pos_{type_name}_R', 'P': f'Pos_{type_name}_P', 'N': f'Pos_{type_name}_N', 'Adjusted_P_Value': f'Pos_{type_name}_FDR'})
                        neg_final = neg.rename(columns={'Name': f'Neg_{type_name}_Name', 'R': f'Neg_{type_name}_R', 'P': f'Neg_{type_name}_P', 'N': f'Neg_{type_name}_N', 'Adjusted_P_Value': f'Neg_{type_name}_FDR'})
                        return pos_final.reset_index(drop=True), neg_final.reset_index(drop=True)

                    pg_pos, pg_neg = get_top_tables(df_g, "Gene")
                    pp_pos, pp_neg = get_top_tables(df_p, "Path")
                    
                    pg_pos = annotate_top_tables(pg_pos, "Gene", "Pos_Gene", gene_to_drugs, pathway_to_genes)
                    pg_neg = annotate_top_tables(pg_neg, "Gene", "Neg_Gene", gene_to_drugs, pathway_to_genes)
                    pp_pos = annotate_top_tables(pp_pos, "Path", "Pos_Path", gene_to_drugs, pathway_to_genes)
                    pp_neg = annotate_top_tables(pp_neg, "Path", "Neg_Path", gene_to_drugs, pathway_to_genes)
                    
                    final_df = pd.concat([pg_pos, pg_neg, pp_pos, pp_neg], axis=1)
                    if not final_df.empty:
                        final_df.to_excel(writer, sheet_name=str(grp_name)[:30].replace(':','_'), index=False)
                        sheet_written = True # Update flag if a sheet is successfully written
                except Exception as e:
                    print(f"❌ Error processing subtype '{grp_name}': {e}. Continuing to the next one.")
                    continue
    except Exception as e: print(f"Failed to save file: {e}")

def run_differential_screening(group1_name, s_group1_idx, group2_name, s_group2_idx, df_gene, df_pathway, p_val_thresh=0.05, save_enabled=True, gene_to_drugs=None, pathway_to_genes=None):
    print(f"\nRunning differential analysis for {group1_name} vs {group2_name} (P < {p_val_thresh})...")
    def _scan_diff(idx1, idx2, df_pool, desc):
        res = []
        for col_name, col_data in tqdm(df_pool.items(), total=df_pool.shape[1], desc=desc, ncols=80, leave=False):
            data1 = col_data.loc[col_data.index.intersection(idx1)].dropna()
            data2 = col_data.loc[col_data.index.intersection(idx2)].dropna()
            if len(data1) < 1 or len(data2) < 1: continue
            
            mean1, mean2 = data1.mean(), data2.mean()
            mean_diff = mean1 - mean2
            
            try: t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=False)
            except: p_val = np.nan
                
            if not pd.isna(p_val):
                res.append({'Name': col_name, 'Mean_Diff': mean_diff, 'P_Value': p_val})
        return pd.DataFrame(res)

    def get_top_tables(df_res, type_name):
        if df_res.empty: return pd.DataFrame(), pd.DataFrame()
        df_sig = df_res[df_res['P_Value'] < p_val_thresh].copy()
        
        # Genes significantly higher in group1 (up-regulated)
        up_df = df_sig[df_sig['Mean_Diff'] > 0].sort_values(by='P_Value', ascending=True)
        up_final = pd.DataFrame({
            f'{group1_name}_{type_name}': up_df['Name'],
            'Diff_Value': up_df['Mean_Diff'],
            'P_Value': up_df['P_Value'],
            'Adjusted_P_Value': up_df.get('Adjusted_P_Value', np.nan)
        })

        # Genes significantly higher in group2 (down-regulated relative to group1)
        down_df = df_sig[df_sig['Mean_Diff'] < 0].sort_values(by='P_Value', ascending=True)
        down_final = pd.DataFrame({
            f'{group2_name}_{type_name}': down_df['Name'],
            'Diff_Value': down_df['Mean_Diff'],
            'P_Value': down_df['P_Value'],
            'Adjusted_P_Value': down_df.get('Adjusted_P_Value', np.nan)
        })
        
        return up_final.reset_index(drop=True), down_final.reset_index(drop=True)

    df_g = _scan_diff(s_group1_idx, s_group2_idx, df_gene, "Scanning Genes")
    df_p = _scan_diff(s_group1_idx, s_group2_idx, df_pathway, "Scanning Paths")

    # [FDR] Combine p-values from genes and pathways for FDR correction
    all_results = pd.concat([df_g, df_p], ignore_index=True)
    if FDR_ENABLED:
        if not all_results.empty:
            # [v82.3 Fix] Handle NaN P-values to avoid length mismatch after FDR correction
            not_na_mask = all_results['P_Value'].notna()
            if not_na_mask.any():
                p_values_no_na = all_results.loc[not_na_mask, 'P_Value']
                fdr_values = multipletests(p_values_no_na, method='fdr_bh')[1]
                all_results['Adjusted_P_Value'] = np.nan
                all_results.loc[not_na_mask, 'Adjusted_P_Value'] = fdr_values
            df_g = all_results.iloc[:len(df_g)].copy()
            df_p = all_results.iloc[len(df_g):].copy()

    pg_up, pg_down = get_top_tables(df_g, "Gene")
    pp_up, pp_down = get_top_tables(df_p, "Path")
    
    pg_up = annotate_top_tables(pg_up, "Gene", f"{group1_name}_Gene", gene_to_drugs, pathway_to_genes)
    pg_down = annotate_top_tables(pg_down, "Gene", f"{group2_name}_Gene", gene_to_drugs, pathway_to_genes)
    pp_up = annotate_top_tables(pp_up, "Path", f"{group1_name}_Path", gene_to_drugs, pathway_to_genes)
    pp_down = annotate_top_tables(pp_down, "Path", f"{group2_name}_Path", gene_to_drugs, pathway_to_genes)
    
    if save_enabled:
        fname = f"DiffScreening_{clean_filename(group1_name)}_vs_{clean_filename(group2_name)}.xlsx"
        try:
            with pd.ExcelWriter(os.path.join(get_desktop_path(), fname), engine='openpyxl') as writer:
                pd.concat([pg_up, pg_down, pp_up, pp_down], axis=1).to_excel(writer, sheet_name='Sig_Diff', index=False)
            print(f"\n🎉 Differential screening complete! Saved to {fname}")
        except Exception as e: print(f"Failed to save file: {e}")
            
    return pg_up, pg_down, pp_up, pp_down

# ================= 8. Main Program =================
def main():
    global IMG_EXT, IMG_DPI, DIST_PLOT_STYLE, FDR_ENABLED
    print("Starting Analysis Tool (v82.1)...")
    gene_file = 'Gene_expression.csv'; pathway_file = 'Signaling_pathway.csv'; clinical_file = 'Clinical_data.csv'; 
    protein_file = 'protein_zscores_TCGA.csv'; pathway_protein_file = 'Signaling_pathway_Protein.csv';
    knockdown_file = 'Knockdown_Dependency_Score.csv'
    has_gene = os.path.exists(gene_file)
    has_path = os.path.exists(pathway_file)
    has_clinical = os.path.exists(clinical_file)
    has_protein = os.path.exists(protein_file)
    has_path_prot = os.path.exists(pathway_protein_file)
    has_knockdown = os.path.exists(knockdown_file)
    
    df_gene = None
    df_pathway = None
    df_clinical = None
    df_protein = None
    df_pathway_protein = None
    df_knockdown, knockdown_subtypes = None, None
    
    need_core_data = has_gene or has_path or has_protein or has_path_prot or has_knockdown
    
    if need_core_data:
        try:
            print("\n[System Initialization - CSV Mode]")
            if has_gene:
                with Spinner(f"Loading {gene_file} ..."): df_gene = load_gene_file(gene_file)
                if df_gene is not None: print(f"   - RNA Gene Expression: {df_gene.shape[0]} patients, {df_gene.shape[1]} genes")
            if has_path:
                with Spinner(f"Loading {pathway_file} ..."): df_pathway = load_pathway_file(pathway_file)
                if df_pathway is not None: print(f"   - RNA Pathway Activity: {df_pathway.shape[0]} patients, {df_pathway.shape[1]} pathways")
            if has_protein:
                with Spinner(f"Loading {protein_file} ..."): df_protein = load_protein_file(protein_file)
                if df_protein is not None: print(f"   - Protein Expression: {df_protein.shape[0]} patients, {df_protein.shape[1]} proteins")
            if has_path_prot:
                with Spinner(f"Loading {pathway_protein_file} ..."): df_pathway_protein = load_pathway_file(pathway_protein_file)
                if df_pathway_protein is not None: print(f"   - Protein Pathway Activity: {df_pathway_protein.shape[0]} patients, {df_pathway_protein.shape[1]} pathways")
            if has_knockdown:
                with Spinner(f"Loading {knockdown_file} ..."):
                    df_knockdown, knockdown_subtypes = load_knockdown_file(knockdown_file)
                if df_knockdown is not None: print(f"   - Gene Knockdown Data: {df_knockdown.shape[0]} cell lines, {df_knockdown.shape[1]} genes")
            if has_clinical:
                with Spinner(f"Loading {clinical_file} ..."): df_clinical = load_clinical_file(clinical_file)
                if df_clinical is not None: print(f"   - Clinical Data: {df_clinical.shape[0]} patients, {df_clinical.shape[1]} features")

            print("\n★ Core system ready! Data will be automatically aligned during analysis.")
        except Exception as e: 
            print(f"\n❌ Core loading error: {e}")
            need_core_data = False
    else:
        print("\n[System Tip] Gene_expression.csv or Pathway file not found. Only [E] External Plotting functions are available.")

    gene_to_drugs = {}
    pathway_to_genes = {}
    if need_core_data:
        print("\n[System Initialization - Loading Annotation Data]")
        gene_to_drugs = build_gene_drug_map('interactions.csv')
        pathway_to_genes = build_pathway_gene_map('Full_Gene_List.csv')
        if gene_to_drugs:
            print(f" -> Successfully loaded drug target info for {len(gene_to_drugs)} genes.")
        if pathway_to_genes:
            print(f" -> Successfully loaded pathway gene compositions for {len(pathway_to_genes)} pathways.")

    subtype_col = None
    if df_clinical is not None and need_core_data:
        print("\n" + "="*60)
        print(" 📡 [Clinical Feature Radar] The system has detected the following columns. You can copy them directly for grouping or analysis:")
        col_names = df_clinical.columns.tolist()
        for i in range(0, len(col_names), 5):
            print(" | ".join(col_names[i:i+5]))
        print("="*60)

        print("\n【Global Stratification Setup】Define multi-level groups. Enter column names (Press Enter to skip/finish):")
        selected_global_cols = []
        while True:
            level = len(selected_global_cols) + 1
            prompt = f"【Stratification {level}: Primary】Enter grouping column (e.g., PAM50), press Enter to skip: " if level == 1 else f"【Stratification {level}: Secondary】Sub-divide within {'_'.join(selected_global_cols)}? Enter column name (Press Enter to skip/finish): "
            c = input(prompt).strip()
            if not c:
                break
            if c in df_clinical.columns:
                selected_global_cols.append(c)
            else:
                print(f"⚠️ Column '{c}' not found. Please re-enter.")
        
        if selected_global_cols:
            if len(selected_global_cols) == 1:
                subtype_col = selected_global_cols[0] # No change needed if only one
                print(f"-> 🎯 Single-level grouping enabled: {subtype_col}")
            else:
                subtype_col = "_".join(selected_global_cols)
                mask = df_clinical[selected_global_cols[0]].notna()
                for col in selected_global_cols[1:]:
                    mask &= df_clinical[col].notna()
                df_clinical.loc[mask, subtype_col] = df_clinical.loc[mask, selected_global_cols[0]].astype(str)
                for col in selected_global_cols[1:]:
                    df_clinical.loc[mask, subtype_col] += "_" + df_clinical.loc[mask, col].astype(str)
                print(f"-> 🎯 {len(selected_global_cols)}-level stratification engine enabled! Grouping by: {subtype_col}")
        else:
            print("-> No global grouping will be applied.")

    print("\n【Global Plot Settings】")
    print("Select the output format for all plots (Press Enter for default: png):")
    print("  (1) png  [Default] (Small file size, good for previews and presentations)")
    print("  (2) pdf  (Vector format, ideal for publications)")
    print("  (3) tiff (Lossless high-resolution bitmap, for journal submissions)")
    print("  (4) svg  (Vector format, for web and post-editing)")
    fmt_choice = input(">> ").strip()
    if fmt_choice == '2': IMG_EXT = "pdf"
    elif fmt_choice == '3': IMG_EXT = "tiff"
    elif fmt_choice == '4': IMG_EXT = "svg"
    else: IMG_EXT = "png"
    
    print("\nEnter plot resolution in DPI (Press Enter for default: 300; high-tier journals often require 600 or 1200):")
    dpi_input = input(">> ").strip()
    try:
        IMG_DPI = int(dpi_input) if dpi_input else 300
    except ValueError:
        IMG_DPI = 300
        
    DIST_PLOT_STYLE = "3"
    style_desc = "Box+Violin fusion plot"
        
    print(f"\n-> 🎯 System locked: Plot format .{IMG_EXT}, Resolution {IMG_DPI} DPI, Style {style_desc}.")
    print(f"-> 📁 [Desktop Protection] Output files will be saved to: Desktop / {SESSION_FOLDER_NAME} /")
    print(f"-> 📈 FDR (False Discovery Rate) correction is enabled by default.")

    while True:
        print("\n" + "="*70)
        print(f"【 Gene Analysis Tool - v82.1 】")
        print(" 💡 [Tip] The system has virtual memory. Signatures you create can be used directly in sections [A], [B], [C]!")
        print(" 💡 [Formula] Supports Ratio (e.g., BCL2/CASP8) or Difference (e.g., BCL2 - CASP8) inputs.")
        
        if need_core_data:
            print("\n[A] Correlation & Grouping Analysis")
            print("  1. Two-Gene Correlation Scatter Plot")
            print("  2. Batch Correlation Matrix (Gene vs Gene Heatmap)")
            print("  3. Patient Co-expression Heatmap")
            print("  4. Single-Target Patient Waterfall Plot")
            print("  5. Target High/Low Group Comparison")

            print("\n[B] Clinical & Survival Analysis")
            print("  6. Clinical Feature Distribution (Box+Violin fusion)")
            print("  7. Single-Gene/Feature Survival Curve (Single KM Plot)")
            print("  8. Dual-Gene/Feature Survival Analysis (Dual-Target KM Plot)")
            print("  9. Hazard Assessment (Cox Forest Plot)")

            print("\n[C] Big Data Screening")
            print("  10. Significant Correlation Screening")
            print("  11. Significant Differential Expression Screening")

            print("\n[D] Advanced Feature Engineering")
            print("  12. Create and Store a Virtual Signature")

        print("\n[E] External Data Plotting (Enrichment, Differential Analysis)")
        print("  13. Single Bar/Lollipop Plot")
        print("  14. Dual Comparison Plot")
        print("  15. Volcano Plot for Differential Analysis")
        print("  16. External Data Waterfall Plot")

        if has_knockdown:
            print("\n[F] CRISPR Knockdown Dependency Analysis")
            print("  17. Dependency Score Waterfall Plot")

        print("\n(Enter 'q' to quit)")
        print("="*70)
        
        mode = input("Enter your choice: ").strip()
        if mode == 'q': break

        try:
            # === [E] External Plotting Block ===
            if mode in ['13', '14', '15', '16']:
                if mode in ['13', '14']:
                    top_n = int(get_default_input("➤ Show Top N items", "15"))
                    chart_style = get_default_input("➤ Style (bar/lollipop)", "bar").lower()
                    color_high = get_default_input("➤ Highlight color (HEX)", "#D62728")
                    x_col = get_default_input("\n➤ Column name for X-axis (default: mlog10p, auto-converted from 'pvalue' column)", "mlog10p")
                    x_label_default = '-log10(P-value)' if x_col == "mlog10p" else x_col
                    x_label = get_default_input("➤ Display name for X-axis", x_label_default)
                    val_col = get_default_input("\n➤ Column name for extra data label (Press Enter to skip)", "")
                    val_label = get_default_input("➤ Prefix for the data label", val_col) if val_col else ""
                    
                    size_col = None
                    size_range = (200, 800)
                    if chart_style == "lollipop":
                        size_col = get_default_input("\n➤ Column name to determine dot size (e.g., NES_abs) [Press Enter to skip]", "")
                        if size_col:
                            s_min = int(get_default_input("   ➤ Minimum size [Default: 200]", "200"))
                            s_max = int(get_default_input("   ➤ Maximum size [Default: 800]", "800"))
                            size_range = (s_min, s_max)

                    if mode == '13':
                        f1 = get_default_input("➤ Data filename (including extension)", "enrichment.csv")
                        title = get_default_input("➤ Chart title (enter 'none' for no title)", "Enrichment Analysis")
                        if title.lower() == 'none': title = ""
                        targets = get_input_list("\nPathways/Genes to highlight? (Highlighted items will be colored)", allow_none=True)
                        df1 = load_ext_enrichment_data(f1, top_n, sort_col=x_col)
                        if df1 is None: continue
                        
                        # [v82.2 Fix] Dynamically adjust figure height and font size to prevent label overlap
                        num_items = len(df1)
                        fig_height = min(40.0, max(6.0, num_items * 0.45))
                        
                        setup_plot_style()
                        fig, ax = plt.subplots(figsize=(10, fig_height))
                        hl_set = set([normalize_name(t) for t in targets])
                        color_norm = "#E0E0E0" if hl_set else "#4C72B0"
                        draw_ext_bar_lollipop(ax, df1, title, hl_set, color_high, color_norm, chart_style, val_col, val_label, x_col, x_label, size_col=size_col, size_range=size_range)
                        out_name = os.path.join(get_desktop_path(), f"Figure_Single_{chart_style}_Top{top_n}.{IMG_EXT}")
                        plt.savefig(out_name, dpi=IMG_DPI, bbox_inches='tight')
                        plt.close()
                        
                        play_beep() # Beep on completion
                        print(f"🎉 Plotting complete! File saved in the session folder.")
                        
                    elif mode == '14':
                        f1 = get_default_input("➤ Filename for Set 1", "enrichment1.csv")
                        title1 = get_default_input("➤ Chart title for Set 1 (enter 'none' for no title)", "Discovery Set")
                        if title1.lower() == 'none': title1 = ""
                        f2 = get_default_input("➤ Filename for Set 2", "enrichment2.csv")
                        title2 = get_default_input("➤ Chart title for Set 2 (enter 'none' for no title)", "Validation Set")
                        if title2.lower() == 'none': title2 = ""
                        main_title = get_default_input("➤ Main chart title (enter 'none' for no title)", "Cross-Dataset Validation")
                        if main_title.lower() == 'none': main_title = ""
                        df1 = load_ext_enrichment_data(f1, top_n, sort_col=x_col)
                        df2 = load_ext_enrichment_data(f2, top_n, sort_col=x_col)
                        if df1 is None or df2 is None: continue

                        # [v82.2 Fix] Dynamically adjust figure height and font size
                        num_items = max(len(df1), len(df2))
                        fig_height = min(40.0, max(6.0, num_items * 0.45))

                        setup_plot_style()
                        overlaps = set(df1['Item_Name'].apply(normalize_name)) & set(df2['Item_Name'].apply(normalize_name))
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, fig_height))
                        draw_ext_bar_lollipop(ax1, df1, title1, overlaps, color_high, "#E0E0E0", chart_style, val_col, val_label, x_col, x_label, size_col=size_col, size_range=size_range)
                        draw_ext_bar_lollipop(ax2, df2, title2, overlaps, color_high, "#E0E0E0", chart_style, val_col, val_label, x_col, x_label, size_col=size_col, size_range=size_range)
                        if main_title:
                            fig.suptitle(f"{main_title}: Top {top_n}", fontsize=22, fontweight='bold', y=0.98)
                        out_name = os.path.join(get_desktop_path(), f"Figure_Dual_{chart_style}_Top{top_n}.{IMG_EXT}")
                        fig.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to prevent title overlap
                        plt.savefig(out_name, dpi=IMG_DPI, bbox_inches='tight')
                        play_beep() # Beep on completion
                        plt.close()
                        print(f"🎉 Plotting complete! File saved in the session folder.")

                elif mode == '15':
                    f1 = get_default_input("➤ Data filename", "volcano.csv")
                    df1 = load_ext_data_basic(f1)
                    if df1 is None: continue
                    
                    x_col = get_default_input("\n➤ Column name for X-axis data (e.g., FoldChange, Correlation)", "Correlation")
                    if x_col not in df1.columns:
                        print(f"[Error] Column not found in CSV: {x_col}"); continue
                    x_trans = get_default_input("➤ X-axis transformation (1: None, 2: log2, 3: -log10)", "1")
                    x_name = get_default_input(f"➤ Display name for X-axis", x_col)

                    y_col = get_default_input("\n➤ Column name for Y-axis data (e.g., pvalue, FDR)", "pvalue")
                    if y_col not in df1.columns:
                        print(f"[Error] Column not found in CSV: {y_col}"); continue
                    y_trans = get_default_input("➤ Y-axis transformation (1: None, 2: log2, 3: -log10) [Recommended for Volcano: 3]", "3")
                    y_name = get_default_input(f"➤ Display name for Y-axis", f"-log10({y_col})" if y_trans == '3' else y_col)
                    
                    s_x = pd.to_numeric(df1[x_col], errors='coerce').astype(float)
                    if x_trans == '2': df1['Plot_X'] = np.log2(s_x.replace(0, np.nan))
                    elif x_trans == '3': 
                        min_pos = s_x[s_x > 0].min() if not s_x[s_x > 0].empty else 1e-50
                        s_x[s_x == 0] = min_pos * 0.1
                        df1['Plot_X'] = -np.log10(s_x)
                    else: df1['Plot_X'] = s_x

                    s_y = pd.to_numeric(df1[y_col], errors='coerce').astype(float)
                    if y_trans == '2': df1['Plot_Y'] = np.log2(s_y.replace(0, np.nan))
                    elif y_trans == '3': 
                        min_pos = s_y[s_y > 0].min() if not s_y[s_y > 0].empty else 1e-50
                        zero_mask = s_y == 0
                        if zero_mask.sum() > 0:
                            s_y[zero_mask] = min_pos * np.random.uniform(0.05, 0.15, size=zero_mask.sum())
                        df1['Plot_Y'] = -np.log10(s_y)
                    else: df1['Plot_Y'] = s_y

                    if x_trans == '2':
                        print("\n[System] log2 transformation detected. Recommended thresholds: 0.585 / -0.585")
                        x_thresh_r = float(get_default_input("➤ X-axis significance threshold [Positive/Up]", "0.585"))
                        x_thresh_l = float(get_default_input("➤ X-axis significance threshold [Negative/Down]", "-0.585"))
                    elif x_trans == '1':
                        if any(k in x_col.lower() for k in ['hr', 'fc', 'ratio', 'fold']):
                            print("\n[System] Ratio data without log2 detected. Recommended thresholds: 1.5 / 0.67")
                            x_thresh_r = float(get_default_input("➤ X-axis significance threshold [Positive/High]", "1.5"))
                            x_thresh_l = float(get_default_input("➤ X-axis significance threshold [Negative/Low]", "0.67"))
                        else:
                            print("\n[System] General numeric data detected. Recommended thresholds: 0.3 / -0.3")
                            x_thresh_r = float(get_default_input("➤ X-axis significance threshold [Positive]", "0.3"))
                            x_thresh_l = float(get_default_input("➤ X-axis significance threshold [Negative]", "-0.3"))
                    else:
                        x_thresh_r = float(get_default_input("➤ X-axis significance threshold [Positive]", "0.5"))
                        x_thresh_l = float(get_default_input("➤ X-axis significance threshold [Negative]", "-0.5"))

                    y_thresh_def = "1.301" if y_trans == '3' else "0.05"
                    y_thresh = float(get_default_input(f"\n➤ Y-axis significance threshold (e.g., -log10(0.05) is ~1.301)", y_thresh_def))
                    y_dir = get_default_input("➤ Y-axis significance direction (1: greater than, 2: less than)", "1" if y_trans == '3' else "2")
                            
                    targets = get_input_list("\nPathways/Genes to label on the volcano plot?", allow_none=True)
                    
                    super_targets = []
                    if targets:
                        super_targets = get_input_list("\n✨【Advanced】From the list above, any to 'super-highlight'? (gold background, red text)", allow_none=True)
                        
                    title = get_default_input("\n➤ Chart title (enter 'none' for no title)", "Volcano Plot")
                    if title.lower() == 'none': title = ""
                    setup_plot_style()
                    
                    n_t = len(targets) if targets else 0
                    fig_w = min(18, 10 + (n_t // 15) * 1.5)
                    fig_h = min(14, 9 + (n_t // 20) * 1.0)
                    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                    
                    if targets: print("⏳ Running collision algorithm to prevent label overlap, please wait...")
                    draw_ext_volcano(ax, df1, title, targets, y_thresh, y_dir, x_thresh_r, x_thresh_l, x_name, y_name, super_targets=super_targets)
                    out_name = os.path.join(get_desktop_path(), f"Figure_Volcano.{IMG_EXT}")
                    plt.tight_layout(rect=[0, 0.03, 1, 0.90])
                    plt.savefig(out_name, dpi=IMG_DPI, bbox_inches='tight')
                    plt.close()
                    play_beep() # Beep on completion
                    print(f"🎉 Plotting complete! File saved in the session folder.")

                elif mode == '16':
                    f1 = get_default_input("➤ Data filename", "enrichment.csv")
                    df1 = load_ext_data_basic(f1) 
                    if df1 is None: continue

                    val_col = get_default_input("➤ Column name for Y-axis values (e.g., NES, logFC, Correlation)", "NES")
                    if val_col not in df1.columns:
                        print(f"[Error] Column not found in CSV: {val_col}"); continue
                    
                    top_n = int(get_default_input("➤ Plot Top N by absolute value (enter 0 for all)", "30"))
                    
                    df1[val_col] = pd.to_numeric(df1[val_col], errors='coerce')
                    df1 = df1.dropna(subset=[val_col])
                    df1['abs_val'] = df1[val_col].abs()
                    
                    if top_n > 0 and top_n < len(df1):
                        df_plot = df1.nlargest(top_n, 'abs_val').sort_values(val_col, ascending=True).reset_index(drop=True)
                    else:
                        df_plot = df1.sort_values(val_col, ascending=True).reset_index(drop=True)
                        top_n = len(df_plot)

                    val_max = df_plot[val_col].max() if not df_plot.empty else 0
                    val_min = df_plot[val_col].min() if not df_plot.empty else 0
                    val_range = val_max - val_min if val_max != val_min else 1.0

                    color_pos = get_default_input("\n➤ Color for positive values (>0)", "#E63946")
                    color_neg = get_default_input("➤ Color for negative values (<0)", "#023E8A")
                    orientation = get_default_input("\n➤ Plot orientation (1: Vertical [Default], 2: Horizontal)", "1")
                    
                    targets = get_input_list("\nPathways/Genes to highlight? (Highlighted items keep their color, others are grayed out)", allow_none=True)
                    hl_set = set([normalize_name(t) for t in targets])

                    print("\n" + "-"*20 + " Adjust Value Axis Range " + "-"*20)
                    print(f"Detected data range is from {val_min:.2f} to {val_max:.2f}.")
                    expand_factor = 1.6 if hl_set else 0.2
                    default_min = float(val_min - val_range * expand_factor)
                    default_max = float(val_max + val_range * expand_factor)
                    
                    user_min_str = input(f"➤ Enter minimum for value axis [Default: {default_min:.2f}]: ").strip()
                    user_min = float(user_min_str) if user_min_str else default_min

                    user_max_str = input(f"➤ Enter maximum for value axis [Default: {default_max:.2f}]: ").strip()
                    user_max = float(user_max_str) if user_max_str else default_max
                    print("-" * (44 + len(" Adjust Value Axis Range ")))
                    
                    title = get_default_input("\n➤ Chart title (enter 'none' for no title)", "Waterfall Plot")
                    if title.lower() == 'none': title = ""
                    setup_plot_style()
                    
                    is_horizontal = (orientation == '2')

                    if is_horizontal:
                        fig_height = min(24, max(8, len(df_plot) * 0.3))
                        fig, ax = plt.subplots(figsize=(10, fig_height))
                    else:
                        fig_width = min(24, max(10, len(df_plot) * 0.2))
                        fig, ax = plt.subplots(figsize=(fig_width, 6))

                    # Call the new drawing function
                    draw_external_waterfall(ax, df_plot, val_col, title, hl_set, color_pos, color_neg, is_horizontal, user_min, user_max)
                    
                    out_name = os.path.join(get_desktop_path(), f"Figure_Waterfall_{val_col}_Top{top_n}.{IMG_EXT}")
                    plt.tight_layout()
                    plt.savefig(out_name, dpi=IMG_DPI, bbox_inches='tight')
                    plt.close()
                    play_beep() # Beep on completion
                    print(f"🎉 Waterfall plot complete! File saved in the session folder.")

            elif mode == '17':
                if df_knockdown is None:
                    print("❌ Error: Knockdown_Dependency_Score.csv not found. This function is unavailable.")
                    continue

                avail_subs = sorted(list(knockdown_subtypes.dropna().unique()))
                print(f"\nAvailable subtypes: {', '.join(avail_subs)}")
                sub_sel = input("➤ Enter subtypes to run (comma-separated, 'all' for all, Enter for All_Samples only): ").strip()
                
                subtypes_to_run = []
                if sub_sel.lower() == 'all':
                    subtypes_to_run = [None] + avail_subs
                elif sub_sel:
                    subtypes_to_run = [s.strip() for s in re.split(r'[,;\t\n|]+', sub_sel) if s.strip() in avail_subs]
                    if not subtypes_to_run: subtypes_to_run = [None]
                else:
                    subtypes_to_run = [None]

                targets = get_input_list("\nGenes to highlight? (Highlighted items will be colored)", allow_none=True)
                hl_set = set([normalize_name(t) for t in targets])

                only_show_hl = False
                if hl_set:
                    print("\nSelect display mode:")
                    print("  (1) Show in context of all genes (keeps background genes for ranking) [Default]")
                    print("  (2) Show only highlighted genes (hides background genes)")
                    ans = input(">> ").strip()
                    if ans == '2': only_show_hl = True

                if not only_show_hl:
                    top_n = int(get_default_input("\n➤ Plot Top N by absolute value (enter 0 for all)", "50"))
                else:
                    top_n = 0

                color_pos = get_default_input("\n➤ Color for positive values (>0)", "#E63946")
                color_neg = get_default_input("➤ Color for negative values (<0)", "#023E8A")
                orientation = get_default_input("\n➤ Plot orientation (1: Vertical [Default], 2: Horizontal)", "1")
                is_horizontal = (orientation == '2')

                for sub in subtypes_to_run:
                    if sub is None:
                        sub_df = df_knockdown
                        grp_name = "All_CellLines"
                    else:
                        cell_lines = knockdown_subtypes[knockdown_subtypes == sub].index
                        sub_df = df_knockdown.loc[cell_lines]
                        grp_name = f"Subtype_{sub}"
                    
                    if sub_df.empty:
                        print(f"Subtype {sub or 'All'} has no data, skipping.")
                        continue

                    print(f"\n---> Processing: {grp_name}")
                    mean_scores = sub_df.mean().sort_values(ascending=True)
                    df_plot = pd.DataFrame({'Item_Name': mean_scores.index, 'Dependency_Score': mean_scores.values})

                    if only_show_hl:
                        hl_mask = df_plot['Item_Name'].apply(lambda x: normalize_name(x) in hl_set)
                        df_plot = df_plot[hl_mask].sort_values('Dependency_Score', ascending=True).reset_index(drop=True)
                        if df_plot.empty:
                            print(f"⚠️ No highlighted genes found in subtype {grp_name}, skipping plot.")
                            continue
                    else:
                        if top_n > 0 and top_n < len(df_plot):
                            df_plot['abs_val'] = df_plot['Dependency_Score'].abs()
                            if hl_set:
                                hl_mask = df_plot['Item_Name'].apply(lambda x: normalize_name(x) in hl_set)
                                df_hl = df_plot[hl_mask]
                                df_others = df_plot[~hl_mask]
                                rem_n = max(0, top_n - len(df_hl))
                                df_top = df_others.nlargest(rem_n, 'abs_val')
                                df_plot = pd.concat([df_hl, df_top]).sort_values('Dependency_Score', ascending=True).reset_index(drop=True)
                            else:
                                df_plot = df_plot.nlargest(top_n, 'abs_val').sort_values('Dependency_Score', ascending=True).reset_index(drop=True)

                    title = f"Gene Dependency Score ({grp_name})"
                    val_col = 'Dependency_Score'
                    
                    setup_plot_style()
                    
                    # Limit max figure size to prevent Memory Error (max width/height 60 inches)
                    if is_horizontal:
                        fig_h = min(60.0, max(8.0, len(df_plot) * 0.25))
                        fig, ax = plt.subplots(figsize=(12, fig_h))
                    else:
                        fig_w = min(60.0, max(10.0, len(df_plot) * 0.15))
                        fig, ax = plt.subplots(figsize=(fig_w, 8))
                        
                    draw_external_waterfall(ax, df_plot, val_col, title, hl_set, color_pos, color_neg, is_horizontal)
                    
                    out_name = os.path.join(get_desktop_path(), f"Figure_Dependency_Waterfall_{clean_filename(grp_name)}.{IMG_EXT}")
                    plt.tight_layout()
                    plt.savefig(out_name, dpi=IMG_DPI, bbox_inches='tight')
                    plt.close()
                    
                    if hl_set:
                        try:
                            df_hl_export = df_plot[df_plot['Item_Name'].apply(lambda x: normalize_name(x) in hl_set)].copy()
                            if 'abs_val' in df_hl_export.columns:
                                df_hl_export = df_hl_export.drop(columns=['abs_val'])
                            if not df_hl_export.empty:
                                excel_name = os.path.join(get_desktop_path(), f"Dependency_Score_Highlighted_{clean_filename(grp_name)}.xlsx")
                                df_hl_export.to_excel(excel_name, index=False)
                        except Exception as e:
                            print(f"⚠️ Failed to export Excel: {e}")
                
                play_beep()
                print(f"🎉 Gene dependency analysis complete! Files saved in the session folder.")

            if not need_core_data and mode.isdigit() and int(mode) < 13:
                print("[System] Core gene library not loaded. Only external plotting modules ([E] section 13-16) are available.")
                continue

            # === [D] Feature Engineering Block ===
            if mode == '12':
                genes_mem = get_input_list("\nEnter gene list for the Signature:")
                if not genes_mem: continue
                l_mode = '1'
                if has_protein or has_path_prot:
                    l_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g = df_protein if l_mode == '2' else df_gene
                df_p = df_pathway_protein if l_mode == '2' else df_pathway
                
                name_mem = input("Name this Signature (e.g., SigA): ").strip() or "SigA"
                if l_mode == '2' and not name_mem.endswith("(Protein)"):
                    name_mem += " (Protein)"
                    
                calc_method = ask_transform_method()
                score_mem, v_mem = calculate_signature_score(genes_mem, df_g, method=calc_method) 
                if score_mem is None: print("❌ Error: No valid genes found"); continue
                df_p[name_mem] = score_mem
                
                try:
                    fname_csv = f"Signature_Data_{clean_filename(name_mem)}.csv"
                    df_out = pd.DataFrame(score_mem, columns=[f"{name_mem}_Score"])
                    if subtype_col: df_out = df_out.join(df_clinical[[subtype_col]], how='left')
                    df_out.to_csv(os.path.join(get_desktop_path(), fname_csv))
                except: pass

                print(f"\n🎉 [Success] Virtual signature '{name_mem}' (integrating {len(v_mem)} genes) has been stored in memory and exported to CSV!")
                print(f"   💡 You can now use '{name_mem}' as a target in options 1-11.")
                play_beep() # Beep on completion
                continue

            # === [A] Correlation Block ===
            elif mode == '1': 
                list1 = get_input_list("\nEnter list 1 (Y-axis, can be multiple):")
                if not list1: continue
                list2 = get_input_list("\nEnter list 2 (X-axis, can be multiple):")
                if not list2: continue
                
                print("\nSelect comparison mode:")
                print("  (1) RNA vs RNA [Default]")
                print("  (2) RNA (Y-axis) vs Protein (X-axis)")
                print("  (3) Protein (Y-axis) vs RNA (X-axis)")
                print("  (4) Protein vs Protein")
                comp_mode = input(">> ").strip() or '1'
                if comp_mode in ['2', '3', '4'] and not has_protein and not has_path_prot:
                    print("⚠️ Warning: Protein data not found. Forcing to RNA vs RNA mode.")
                    comp_mode = '1'

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                
                if comp_mode == '2':
                    df_g1, df_p1 = df_gene, df_pathway
                    df_g2, df_p2 = df_protein, df_pathway_protein
                    suffix1, suffix2 = " (RNA)", " (Protein)"
                elif comp_mode == '3':
                    df_g1, df_p1 = df_protein, df_pathway_protein
                    df_g2, df_p2 = df_gene, df_pathway
                    suffix1, suffix2 = " (Protein)", " (RNA)"
                elif comp_mode == '4':
                    df_g1, df_p1 = df_protein, df_pathway_protein
                    df_g2, df_p2 = df_protein, df_pathway_protein
                    suffix1, suffix2 = " (Protein)", " (Protein)"
                else:
                    df_g1, df_p1 = df_gene, df_pathway
                    df_g2, df_p2 = df_gene, df_pathway
                    suffix1, suffix2 = "", ""

                for n1 in list1:
                    for n2 in list2:
                        d1, t1 = get_data_by_name(n1, df_g1, df_p1)
                        d2, t2 = get_data_by_name(n2, df_g2, df_p2)
                        if d1 is None or d2 is None: continue
                        d1 = transform_series(d1, calc_method)
                        d2 = transform_series(d2, calc_method)
                        
                        display_n1 = n1 + suffix1 if suffix1 and 'Protein' not in n1 else n1
                        display_n2 = n2 + suffix2 if suffix2 else n2
                        run_subtype_correlation(display_n1, d1, display_n2, d2, df_clinical, local_sub_col, plot_enabled=True, subtypes_to_run=subs_to_run)
                print("\n🎉 Batch scatter plot generation complete!")
                play_beep() # Beep on completion

            elif mode == '2': 
                list1 = get_input_list("\nEnter list 1 (Y-axis):")
                if not list1: continue
                list2 = get_input_list("\nEnter list 2 (X-axis):")
                if not list2: continue
                
                print("\nSelect comparison mode:")
                print("  (1) RNA vs RNA [Default]")
                print("  (2) RNA (Y-axis) vs Protein (X-axis)")
                print("  (3) Protein (Y-axis) vs RNA (X-axis)")
                print("  (4) Protein vs Protein")
                comp_mode = input(">> ").strip() or '1'
                if comp_mode in ['2', '3', '4'] and not has_protein and not has_path_prot:
                    print("⚠️ Warning: Protein data not found. Forcing to RNA vs RNA mode.")
                    comp_mode = '1'
                    
                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                cmap = ask_for_palette('diverging')
                
                if comp_mode == '2':
                    df_g1, df_p1 = df_gene, df_pathway
                    df_g2, df_p2 = df_protein, df_pathway_protein
                    suffix1, suffix2 = " (RNA)", " (Protein)"
                elif comp_mode == '3':
                    df_g1, df_p1 = df_protein, df_pathway_protein
                    df_g2, df_p2 = df_gene, df_pathway
                    suffix1, suffix2 = " (Protein)", " (RNA)"
                elif comp_mode == '4':
                    df_g1, df_p1 = df_protein, df_pathway_protein
                    df_g2, df_p2 = df_protein, df_pathway_protein
                    suffix1, suffix2 = " (Protein)", " (Protein)"
                else:
                    df_g1, df_p1 = df_gene, df_pathway
                    df_g2, df_p2 = df_gene, df_pathway
                    suffix1, suffix2 = "", ""
                    
                t_g1 = transform_df(df_g1, calc_method)
                t_p1 = transform_df(df_p1, calc_method)
                t_g2 = transform_df(df_g2, calc_method)
                t_p2 = transform_df(df_p2, calc_method)
                
                run_batch_correlation(list1, list2, t_g1, t_p1, t_g2, t_p2, df_clinical, local_sub_col, suffix1, suffix2, subtypes_to_run=subs_to_run, cmap=cmap)
                print("\n🎉 Batch correlation matrix generation complete!")
                play_beep() # Beep on completion
                
            elif mode == '3': 
                t_name = input("\nEnter the base target: ").strip()
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for the base target: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g1 = df_protein if t_mode == '2' else df_gene
                df_p1 = df_pathway_protein if t_mode == '2' else df_pathway
                t_data, t_type = get_data_by_name(t_name, df_g1, df_p1)
                if t_data is None: continue
                disp_t = t_name + " (Protein)" if t_mode == '2' else t_name
                
                list_genes = get_input_list("\nEnter the list of genes/pathways to compare against:")
                if not list_genes: continue
                l_mode = '1'
                if has_protein or has_path_prot:
                    l_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g2 = df_protein if l_mode == '2' else df_gene
                df_p2 = df_pathway_protein if l_mode == '2' else df_pathway
                suffix2 = " (Protein)" if l_mode == '2' else ""
                
                hl_genes = get_input_list("\nGenes/Pathways to highlight? (Highlighted items will be colored)", allow_none=True)
                hl_set = set([normalize_name(g) for g in hl_genes]) if hl_genes else None

                cmap = ask_for_palette('diverging')

                print("\nSelect plot orientation: (1) Horizontal [Default]  (2) Vertical")
                ori = input(">> ").strip() or '1'
                print("\nSelect sort order for comparison genes: (1) Descending [Default]  (2) Ascending")
                sort_ord = input(">> ").strip() or '1'
                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                t_data = transform_series(t_data, calc_method)
                t_gene_2 = transform_df(df_g2, calc_method)
                t_path_2 = transform_df(df_p2, calc_method)
                run_patient_heatmap(disp_t, t_data, list_genes, t_gene_2, t_path_2, df_clinical, local_sub_col, ori, sort_ord, subtypes_to_run=subs_to_run, suffix=suffix2, hl_set=hl_set, cmap=cmap)
                print("\n🎉 Patient co-expression heatmap generation complete!")
                play_beep() # Beep on completion

            elif mode == '4': 
                t_name = input("\nEnter the base target: ").strip()
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for this target: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g = df_protein if t_mode == '2' else df_gene
                df_p = df_pathway_protein if t_mode == '2' else df_pathway
                
                t_data, _ = get_data_by_name(t_name, df_g, df_p)
                if t_data is None: continue
                disp_t = t_name + " (Protein)" if t_mode == '2' else t_name
                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                palette = ask_for_palette('categorical')
                t_data = transform_series(t_data, calc_method)
                run_patient_waterfall(disp_t, t_data, df_clinical, local_sub_col, subtypes_to_run=subs_to_run, palette=palette)
                print("\n🎉 Patient waterfall plot generation complete!")
                play_beep() # Beep on completion

            elif mode == '5': 
                grp_name = input("\nEnter the grouping target: ").strip()
                grp_mode = '1'
                if has_protein or has_path_prot:
                    grp_mode = input(">> Data type for the grouping target: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g1 = df_protein if grp_mode == '2' else df_gene
                df_p1 = df_pathway_protein if grp_mode == '2' else df_pathway
                
                grp_data, _ = get_data_by_name(grp_name, df_g1, df_p1)
                if grp_data is None: continue
                disp_grp = grp_name + " (Protein)" if grp_mode == '2' else grp_name
                
                list_genes = get_input_list("\nEnter the list of genes/pathways to compare:")
                if not list_genes: continue
                l_mode = '1'
                if has_protein or has_path_prot:
                    l_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g2 = df_protein if l_mode == '2' else df_gene
                df_p2 = df_pathway_protein if l_mode == '2' else df_pathway
                
                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                grp_data = transform_series(grp_data, calc_method)
                t_gene2 = transform_df(df_g2, calc_method)
                t_path2 = transform_df(df_p2, calc_method)
                compare_dict = {}
                for g in list_genes:
                    d, _ = get_data_by_name(g, t_gene2, t_path2)
                    if d is not None: 
                        disp_g = g + " (Protein)" if l_mode == '2' else g
                        compare_dict[disp_g] = d
                if compare_dict:
                    run_high_low_comparison(disp_grp, grp_data, compare_dict, df_clinical, local_sub_col, subtypes_to_run=subs_to_run)
                    print("\n🎉 High/Low group comparison complete!")
                    play_beep() # Beep on completion

            # === [B] Clinical & Survival Block ===
            elif mode == '6':
                if not has_clinical:
                    print("❌ Error: Clinical_data.csv not found. This function is unavailable.")
                    continue
                targets = get_input_list("Enter target names:")
                if not targets: continue
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g = df_protein if t_mode == '2' else df_gene
                df_p = df_pathway_protein if t_mode == '2' else df_pathway

                print("\nAvailable clinical features:")
                col_names = df_clinical.columns.tolist()
                for i in range(0, len(col_names), 5):
                    print(" | ".join(col_names[i:i+5]))
                selected_cols = get_input_list("\nEnter primary clinical features for analysis (X-axis) (Type N to skip, or paste from Excel/above):", allow_none=True)

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)

                calc_method = ask_transform_method()
                palette = ask_for_palette('categorical')
                for t in targets:
                    td, dtype = get_data_by_name(t, df_g, df_p)
                    if td is None: continue
                    disp_t = t + " (Protein)" if t_mode == '2' else t
                    td = transform_series(td, calc_method)
                    run_clinical_analysis(disp_t, td, df_clinical, local_sub_col, plot_enabled=True, selected_features=selected_cols, hue_col=None, subtypes_to_run=subs_to_run, palette=palette)
                print("\n🎉 Clinical feature distribution analysis complete!")
                play_beep() # Beep on completion

            elif mode == '7':
                if not has_clinical: continue
                targets = get_input_list("Enter target names:")
                if not targets: continue
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g = df_protein if t_mode == '2' else df_gene
                df_p = df_pathway_protein if t_mode == '2' else df_pathway
                
                survival_pairs = detect_survival_pairs(df_clinical)
                
                selected_pairs = []
                if survival_pairs:
                    print("\nDetected the following survival analysis pairs:")
                    for i, (t_col, s_col) in enumerate(survival_pairs):
                        print(f"  ({i+1}) Time: {t_col}, Status: {s_col}")
                    print(f"  ({len(survival_pairs)+1}) Run all automatically")
                    print(f"  ({len(survival_pairs)+2}) Enter manually")
                    
                    pair_choice = input("Select a pair to analyze (enter number): ").strip()
                    
                    if pair_choice.isdigit():
                        choice_idx = int(pair_choice)
                        if 1 <= choice_idx <= len(survival_pairs):
                            selected_pairs.append(survival_pairs[choice_idx-1])
                        elif choice_idx == len(survival_pairs) + 1:
                            selected_pairs = survival_pairs
                        elif choice_idx == len(survival_pairs) + 2:
                            # Fallback to manual input
                            t_col = input("Time column: ").strip() or 'OS_Months'
                            s_col = input("Status column: ").strip() or 'OS_Status'
                            selected_pairs.append((t_col, s_col))
                        else:
                            print("Invalid choice. Using default (OS_Months, OS_Status).")
                            selected_pairs.append(('OS_Months', 'OS_Status'))
                    else:
                        print("Invalid choice. Using default (OS_Months, OS_Status).")
                        selected_pairs.append(('OS_Months', 'OS_Status'))
                else:
                    print("\nNo standard survival pairs detected. Please enter manually.")
                    t_col = input("Time column: ").strip() or 'OS_Months'
                    s_col = input("Status column: ").strip() or 'OS_Status'
                    selected_pairs.append((t_col, s_col))

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                for t_col, s_col in selected_pairs:
                    for t in targets:
                        td, dtype = get_data_by_name(t, df_g, df_p)
                        if td is None: continue
                        disp_t = t + " (Protein)" if t_mode == '2' else t
                        td = transform_series(td, calc_method)
                        run_survival_analysis(disp_t, td, df_clinical, local_sub_col, t_col, s_col, subtypes_to_run=subs_to_run)
                print("\n🎉 Single-gene survival curve analysis complete!")
                play_beep() # Beep on completion

            elif mode == '8':
                if not has_clinical: continue
                t1_name = input("\nEnter primary base target (Target 1): ").strip()
                t1_mode = '1'
                if has_protein or has_path_prot:
                    t1_mode = input(">> Data type for Target 1: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g1 = df_protein if t1_mode == '2' else df_gene
                df_p1 = df_pathway_protein if t1_mode == '2' else df_pathway
                t1_data, _ = get_data_by_name(t1_name, df_g1, df_p1)
                if t1_data is None: 
                    print(f"Target {t1_name} not found.")
                    continue
                disp_t1 = t1_name + " (Protein)" if t1_mode == '2' else t1_name

                list2 = get_input_list("\nEnter list for Target 2 (supports batch processing):")
                if not list2: continue
                l2_mode = '1'
                if has_protein or has_path_prot:
                    l2_mode = input(">> Data type for Target 2 list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g2 = df_protein if l2_mode == '2' else df_gene
                df_p2 = df_pathway_protein if l2_mode == '2' else df_pathway
                
                survival_pairs = detect_survival_pairs(df_clinical)
                
                selected_pairs = []
                if survival_pairs:
                    print("\nDetected the following survival analysis pairs:")
                    for i, (t_col, s_col) in enumerate(survival_pairs):
                        print(f"  ({i+1}) Time: {t_col}, Status: {s_col}")
                    print(f"  ({len(survival_pairs)+1}) Run all automatically")
                    print(f"  ({len(survival_pairs)+2}) Enter manually")
                    
                    pair_choice = input("Select a pair to analyze (enter number): ").strip()
                    
                    if pair_choice.isdigit():
                        choice_idx = int(pair_choice)
                        if 1 <= choice_idx <= len(survival_pairs):
                            selected_pairs.append(survival_pairs[choice_idx-1])
                        elif choice_idx == len(survival_pairs) + 1:
                            selected_pairs = survival_pairs
                        elif choice_idx == len(survival_pairs) + 2:
                            # Fallback to manual input
                            t_col = input("Time column: ").strip() or 'OS_Months'
                            s_col = input("Status column: ").strip() or 'OS_Status'
                            selected_pairs.append((t_col, s_col))
                        else:
                            print("Invalid choice. Using default (OS_Months, OS_Status).")
                            selected_pairs.append(('OS_Months', 'OS_Status'))
                    else:
                        print("Invalid choice. Using default (OS_Months, OS_Status).")
                        selected_pairs.append(('OS_Months', 'OS_Status'))
                else:
                    print("\nNo standard survival pairs detected. Please enter manually.")
                    t_col = input("Time column: ").strip() or 'OS_Months'
                    s_col = input("Status column: ").strip() or 'OS_Status'
                    selected_pairs.append((t_col, s_col))

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()

                t1_data = transform_series(t1_data, calc_method)
                t_gene_2 = transform_df(df_g2, calc_method)
                t_path_2 = transform_df(df_p2, calc_method)

                print(f"\nStarting dual-gene survival analysis engine (for {len(list2)} pairs)...")
                for t_col, s_col in selected_pairs:
                    for t2_name in list2:
                        t2_data, _ = get_data_by_name(t2_name, t_gene_2, t_path_2)
                        if t2_data is not None:
                            disp_t2 = t2_name + " (Protein)" if l2_mode == '2' else t2_name
                            run_dual_survival_analysis(disp_t1, t1_data, disp_t2, t2_data, df_clinical, local_sub_col, t_col, s_col, subtypes_to_run=subs_to_run)
                print("\n🎉 Batch dual-gene survival analysis complete!")
                play_beep() # Beep on completion

            elif mode == '9': 
                if not has_clinical: continue
                targets = get_input_list("Enter target names:")
                if not targets: continue
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for this list: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g = df_protein if t_mode == '2' else df_gene
                df_p = df_pathway_protein if t_mode == '2' else df_pathway
                
                survival_pairs = detect_survival_pairs(df_clinical)
                
                selected_pairs = []
                if survival_pairs:
                    print("\nDetected the following survival analysis pairs:")
                    for i, (t_col, s_col) in enumerate(survival_pairs):
                        print(f"  ({i+1}) Time: {t_col}, Status: {s_col}")
                    print(f"  ({len(survival_pairs)+1}) Run all automatically")
                    print(f"  ({len(survival_pairs)+2}) Enter manually")
                    
                    pair_choice = input("Select a pair to analyze (enter number): ").strip()
                    
                    if pair_choice.isdigit():
                        choice_idx = int(pair_choice)
                        if 1 <= choice_idx <= len(survival_pairs):
                            selected_pairs.append(survival_pairs[choice_idx-1])
                        elif choice_idx == len(survival_pairs) + 1:
                            selected_pairs = survival_pairs
                        elif choice_idx == len(survival_pairs) + 2:
                            t_col = input("Time column: ").strip() or 'OS_Months'
                            s_col = input("Status column: ").strip() or 'OS_Status'
                            selected_pairs.append((t_col, s_col))
                        else:
                            print("Invalid choice. Using default (OS_Months, OS_Status).")
                            selected_pairs.append(('OS_Months', 'OS_Status'))
                    else:
                        print("Invalid choice. Using default (OS_Months, OS_Status).")
                        selected_pairs.append(('OS_Months', 'OS_Status'))
                else:
                    print("\nNo standard survival pairs detected. Please enter manually.")
                    t_col = input("Time column: ").strip() or 'OS_Months'
                    s_col = input("Status column: ").strip() or 'OS_Status'
                    selected_pairs.append((t_col, s_col))

                layout_choice = '1'
                if len(targets) > 1:
                    print("\nSelect summary plot orientation: (1) Vertical [Default]  (2) Horizontal")
                    layout_choice = input(">> ").strip() or '1'
                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                
                for t_col, s_col in selected_pairs:
                    print(f"\n---> Processing survival pair: {t_col} / {s_col}")
                    batch_results = []
                    all_summary_data = [] 
                    for t in targets:
                        td, dtype = get_data_by_name(t, df_g, df_p)
                        if td is None: continue
                        disp_t = t + " (Protein)" if t_mode == '2' else t
                        td = transform_series(td, calc_method)
                        res = run_cox_analysis(disp_t, td, df_clinical, local_sub_col, t_col, s_col, plot_individual=True, save_excel=False, subtypes_to_run=subs_to_run)
                        if res: 
                            batch_results.append(res)
                            all_summary_data.extend(res)
                    if all_summary_data:
                        summary_df = pd.DataFrame(all_summary_data)
                        cols = ['Gene/Signature', 'Scope', 'Hazard_Ratio', 'CI_Lower_95', 'CI_Upper_95', 'P_Value', 'N', 'Events', 'Time_Col', 'Status_Col']
                        cols = [c for c in cols if c in summary_df.columns]
                        summary_df = summary_df[cols]
                        summary_df.to_excel(os.path.join(get_desktop_path(), f"Batch_Cox_Summary_{clean_filename(t_col)}_{clean_filename(s_col)}.xlsx"), index=False)
                    if len(batch_results) > 1:
                        all_scopes = set()
                        for res_list in batch_results:
                            for item in res_list: all_scopes.add(item['Scope'])
                        for scope in all_scopes:
                            draw_batch_forest_plot(batch_results, scope, orientation=layout_choice, time_col=t_col, status_col=s_col)
                            
                print("\n🎉 Cox Hazard Forest Plot analysis complete!")
                play_beep() # Beep on completion

            # === [C] Big Data Screening Block ===
            elif mode == '10':
                targets = get_input_list("\nEnter base target(s) for screening (multiple allowed):")
                if not targets: continue
                t_mode = '1'
                if has_protein or has_path_prot:
                    t_mode = input(">> Data type for base target(s): (1) RNA [Default]  (2) Protein: ").strip() or '1'
                df_g1 = df_protein if t_mode == '2' else df_gene
                df_p1 = df_pathway_protein if t_mode == '2' else df_pathway
                
                print("\nSelect database to screen against:")
                print("  (1) RNA database [Default]")
                print("  (2) Protein database")
                db_mode = input(">> ").strip() or '1'
                if db_mode == '2' and not has_protein and not has_path_prot:
                    print("⚠️ Protein data not found. Forcing to RNA database.")
                    db_mode = '1'
                df_g2 = df_protein if db_mode == '2' else df_gene
                df_p2 = df_pathway_protein if db_mode == '2' else df_pathway
                
                try:
                    p_val_input = input("➤ Enter P-value threshold [Default: 0.05]: ").strip()
                    p_val_thresh = float(p_val_input) if p_val_input else 0.05
                except ValueError:
                    p_val_thresh = 0.05

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                calc_method = ask_transform_method()
                t_gene_2 = transform_df(df_g2, calc_method) # Screening DB
                t_path_2 = transform_df(df_p2, calc_method) # Screening DB
                
                for t in targets:
                    td, dtype = get_data_by_name(t, df_g1, df_p1) # Target
                    if td is None: continue
                    disp_t = t + " (Protein)" if t_mode == '2' else t
                    td = transform_series(td, calc_method)
                    
                    run_subtype_screening(disp_t, td, t_gene_2, t_path_2, df_clinical, local_sub_col, p_val_thresh, subtypes_to_run=subs_to_run, gene_to_drugs=gene_to_drugs, pathway_to_genes=pathway_to_genes)
                print("\n🎉 Batch correlation screening complete!")
                play_beep() # Beep on completion

            elif mode == '11':
                print("\nSelect grouping basis:")
                print("  (1) By high/low expression of a target gene/pathway")
                if has_clinical: print("  (2) By a clinical feature (e.g., Tumor vs Normal)")
                
                diff_choice = input(">> ").strip()

                try:
                    p_val_input = input("➤ Enter P-value threshold [Default: 0.05]: ").strip()
                    p_val_thresh = float(p_val_input) if p_val_input else 0.05
                except ValueError:
                    p_val_thresh = 0.05

                local_sub_col, subs_to_run = ask_local_subtyping(df_clinical, subtype_col)
                if subs_to_run is None:
                    if local_sub_col and local_sub_col in df_clinical.columns:
                        subs_to_run = [None] + list(df_clinical[local_sub_col].dropna().unique())
                    else:
                        subs_to_run = [None]

                if diff_choice == '1':
                    grp_name = input("\nEnter the grouping target: ").strip()
                    grp_mode = '1'
                    if has_protein or has_path_prot:
                        grp_mode = input(">> Data type for the grouping target: (1) RNA [Default]  (2) Protein: ").strip() or '1'
                    df_g1 = df_protein if grp_mode == '2' else df_gene
                    df_p1 = df_pathway_protein if grp_mode == '2' else df_pathway
                    grp_data, _ = get_data_by_name(grp_name, df_g1, df_p1)
                    if grp_data is None: continue
                    disp_grp = grp_name + " (Protein)" if grp_mode == '2' else grp_name
                    
                    print("\nSelect database to screen against:")
                    print("  (1) RNA database [Default]")
                    print("  (2) Protein database")
                    db_mode = input(">> ").strip() or '1'
                    if db_mode == '2' and not has_protein and not has_path_prot:
                        print("⚠️ Protein data not found. Forcing to RNA database.")
                        db_mode = '1'
                    df_g2 = df_protein if db_mode == '2' else df_gene
                    df_p2 = df_pathway_protein if db_mode == '2' else df_pathway

                    print("Select splitting method:")
                    print("  (1) Median (Top 50% vs Bottom 50%) [Default]")
                    print("  (2) Tertiles (Top 33% vs Bottom 33%)")
                    print("  (3) Quartiles (Top 25% vs Bottom 25%)")
                    q_choice = input(">> ").strip()
                    
                    if q_choice == '2': val_high, val_low = grp_data.quantile(0.6667), grp_data.quantile(0.3333)
                    elif q_choice == '3': val_high, val_low = grp_data.quantile(0.75), grp_data.quantile(0.25)
                    else: val_high = val_low = grp_data.median()
                    
                    idx_high, idx_low = grp_data[grp_data >= val_high].index, (grp_data[grp_data <= val_low].index if val_high > val_low else grp_data[grp_data < val_low].index)
                    
                    for sub in subs_to_run:
                        cur_idx_high, cur_idx_low = idx_high, idx_low
                        name_high, name_low = f"{disp_grp}_High", f"{disp_grp}_Low"
                        if sub is not None:
                            print(f"\n---> Processing subtype: {sub}")
                            idx_sub = df_clinical[df_clinical[local_sub_col] == sub].index
                            cur_idx_high = idx_high.intersection(idx_sub)
                            cur_idx_low = idx_low.intersection(idx_sub)
                            name_high, name_low = f"{disp_grp}_High_in_{sub}", f"{disp_grp}_Low_in_{sub}"
                        
                        run_differential_screening(name_high, cur_idx_high, name_low, cur_idx_low, df_g2, df_p2, p_val_thresh, gene_to_drugs=gene_to_drugs, pathway_to_genes=pathway_to_genes)
                            
                    play_beep() # Beep on completion

                elif diff_choice == '2' and has_clinical:
                    feat = input(f"\nEnter clinical feature column name (Available: {', '.join(df_clinical.columns[:10])}...): ").strip()
                    if feat not in df_clinical.columns: continue
                    cats = df_clinical[feat].dropna().unique()
                    print(f"Detected groups: {', '.join(map(str, cats))}")
                    c1 = input("Enter name for Group 1: ").strip()
                    c2 = input("Enter name for Group 2: ").strip()
                    idx1, idx2 = df_clinical[df_clinical[feat] == c1].index, df_clinical[df_clinical[feat] == c2].index
                    
                    print("\nSelect database to screen against:")
                    print("  (1) RNA database [Default]")
                    print("  (2) Protein database")
                    db_mode = input(">> ").strip() or '1'
                    if db_mode == '2' and not has_protein and not has_path_prot:
                        print("⚠️ Protein data not found. Forcing to RNA database.")
                        db_mode = '1'
                    df_g2 = df_protein if db_mode == '2' else df_gene
                    df_p2 = df_pathway_protein if db_mode == '2' else df_pathway

                    for sub in subs_to_run:
                        cur_idx1, cur_idx2 = idx1, idx2
                        name1, name2 = c1, c2
                        if sub is not None:
                            print(f"\n---> Processing subtype: {sub}")
                            idx_sub = df_clinical[df_clinical[local_sub_col] == sub].index
                            cur_idx1 = idx1.intersection(idx_sub)
                            cur_idx2 = idx2.intersection(idx_sub)
                            name1, name2 = f"{c1}_in_{sub}", f"{c2}_in_{sub}"

                        run_differential_screening(name1, cur_idx1, name2, cur_idx2, df_g2, df_p2, p_val_thresh, gene_to_drugs=gene_to_drugs, pathway_to_genes=pathway_to_genes)
                    play_beep() # Beep on completion

        except Exception as e: print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
