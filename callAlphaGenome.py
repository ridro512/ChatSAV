#!/usr/bin/env python3
"""
callAlphaGenome.py - File for calling AlphaGenome predictions

Input: Variant coordinates in chr:pos:ref:alt format and genome build
Output: AlphaGenome predictions for splice sites, gene expression, and other genomic features

Note: Requires the alphagenome package to be installed:
pip install alphagenome

"""

import os
import pandas as pd
from typing import Dict, Any, Optional, List, Tuple
import warnings

API_KEY = os.getenv("ALPHA_GENOME_API_KEY")
if not API_KEY:
    raise ValueError(
            "AlphaGenome API key not found. Please set the ALPHA_GENOME_API_KEY environment variable.\n"
            "You can do this by:\n"
            "1. Running: export ALPHA_GENOME_API_KEY='your_api_key_here'\n"
            "2. Creating a .env file with: ALPHA_GENOME_API_KEY=your_api_key_here\n"
            "3. Setting it in your shell profile (.bashrc, .zshrc, etc.)"
        )

try:
    from alphagenome import colab_utils
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    ALPHAGENOME_AVAILABLE = True
except ImportError as e:
    print(f"Warning: AlphaGenome package not available: {e}")
    print("To install: pip install alphagenome")
    ALPHAGENOME_AVAILABLE = False

def parse_variant_coordinates(variant_coord: str, hg: str) -> Tuple[str, int, str, str]:
    """
    Parse variant coordinates from string format.
    """
    try:
        chrom, pos, ref, alt = variant_coord.split(":")
        return chrom, int(pos), ref, alt
    except ValueError as e:
        raise ValueError(f"Invalid variant format. Expected 'chr:pos:ref:alt', got '{variant_coord}'") from e

def call_alphagenome(variant_coord: str, hg: str, sequence_length: str = "100KB") -> Dict[str, Any]:
    """
    Query AlphaGenome for comprehensive variant effect predictions.
    """
    results = {
        "variant": variant_coord,
        "genome_build": f"GRCh{hg}",
        "sequence_length": sequence_length,
        "alphagenome_scores": {},
        "api_endpoint": "alphagenome",
        "model_info": {
            "name": "AlphaGenome",
            "provider": "Google DeepMind",
            "version": "2025"
        }
    }
    
    if not ALPHAGENOME_AVAILABLE:
        results["error"] = "AlphaGenome package not installed. Run: pip install alphagenome"
        return results
    
    # currently AlphaGenome supports GRCh38/hg38 primarily
    if hg == "37":
        print("Warning: AlphaGenome is optimized for GRCh38. Results for GRCh37 may be less accurate.")
        organism = dna_client.Organism.HOMO_SAPIENS
    elif hg == "38":
        organism = dna_client.Organism.HOMO_SAPIENS
    else:
        results["error"] = f"Unsupported genome build: GRCh{hg}. AlphaGenome supports GRCh38."
        return results
    
    try:
        chrom, pos, ref, alt = parse_variant_coordinates(variant_coord, hg)
        
        # initialize the AlphaGenome model
        try:
            dna_model = dna_client.create(API_KEY)
        except Exception as e:
            results["error"] = f"Failed to initialize AlphaGenome client: {str(e)}"
            return results
        
        # create variant object
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
        )
        
        # set sequence length
        length_map = {
            "2KB": dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_2KB'],
            "16KB": dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_16KB'], 
            "100KB": dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_100KB'],
            "500KB": dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_500KB'],
            "1MB": dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
        }
        
        if sequence_length not in length_map:
            sequence_length = "100KB"
            
        seq_length = length_map[sequence_length]
        
        # the input interval is derived from the variant
        interval = variant.reference_interval.resize(seq_length)
        
        # score variant using recommended scorers for splice-relevant modalities
        # call each scorer individually to avoid the hashable Variant issue
        all_variant_scores = []
        
        individual_scorers = [
            ("SPLICE_SITES", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('SPLICE_SITES')),
            ("SPLICE_SITE_USAGE", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('SPLICE_SITE_USAGE')),
            ("RNA_SEQ", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('RNA_SEQ')),
            ("CAGE", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('CAGE')),
            ("ATAC", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('ATAC')),
            ("DNASE", variant_scorers.RECOMMENDED_VARIANT_SCORERS.get('DNASE')),
        ]
        
        successful_scorers = []
        
        print(f"Scoring variant {variant_coord} with individual AlphaGenome scorers...")
        
        for scorer_name, scorer in individual_scorers:
            if scorer is None:
                continue
                
            try:
                print(f"  Trying {scorer_name} scorer...")
                variant_scores = dna_model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=[scorer]
                )
                all_variant_scores.append(variant_scores)
                successful_scorers.append(scorer_name)
                print(f"  ✓ {scorer_name} scorer successful")
                
            except Exception as scorer_error:
                print(f"  ✗ {scorer_name} scorer failed: {scorer_error}")
                continue
        
        if not all_variant_scores:
            raise Exception("All individual scorers failed")
        
        print(f"✓ Successfully scored with {len(successful_scorers)} scorers: {', '.join(successful_scorers)}")
        
        # combine all results with better error handling
        combined_df_scores = []
        for i, (scorer_name, variant_score_list) in enumerate(zip(successful_scorers, all_variant_scores)):
            try:
                print(f"  Processing results from {scorer_name}...")
                df_score = variant_scorers.tidy_scores(variant_score_list)
                combined_df_scores.append(df_score)
                print(f"  ✓ {scorer_name} results processed: {len(df_score)} predictions")
            except Exception as process_error:
                print(f"  ✗ Error processing {scorer_name} results: {process_error}")
                continue
        
        if not combined_df_scores:
            raise Exception("Failed to process results from any scorer")
        
        # concatenate all dataframes
        import pandas as pd
        df_scores = pd.concat(combined_df_scores, ignore_index=True)
        
        # process results by output type with error handling
        try:
            print("  Processing splice predictions...")
            splice_results = process_splice_predictions(df_scores)
        except Exception as e:
            print(f"  Warning: Error processing splice predictions: {e}")
            splice_results = {"error": str(e)}
        
        try:
            print("  Processing expression predictions...")
            expression_results = process_expression_predictions(df_scores)
        except Exception as e:
            print(f"  Warning: Error processing expression predictions: {e}")
            expression_results = {"error": str(e)}
        
        try:
            print("  Processing chromatin predictions...")
            chromatin_results = process_chromatin_predictions(df_scores)
        except Exception as e:
            print(f"  Warning: Error processing chromatin predictions: {e}")
            chromatin_results = {"error": str(e)}
        
        try:
            print("  Getting summary statistics...")
            summary_stats = get_summary_statistics(df_scores)
        except Exception as e:
            print(f"  Warning: Error getting summary statistics: {e}")
            summary_stats = {"error": str(e)}
        
        results["alphagenome_scores"] = {
            "splice_predictions": splice_results,
            "expression_predictions": expression_results,
            "chromatin_predictions": chromatin_results,
            "summary_stats": summary_stats,
            "raw_scores_count": len(df_scores)
        }
        
        # add metadata about successful scorers
        results["successful_scorers"] = successful_scorers
        results["total_scorers_attempted"] = len([s for s in individual_scorers if s[1] is not None])
        
        # add top predictions summary with error handling
        try:
            print("  Getting top predictions...")
            results["top_predictions"] = get_top_predictions(df_scores)
        except Exception as e:
            print(f"  Warning: Error getting top predictions: {e}")
            results["top_predictions"] = []
        
        print(f"✓ AlphaGenome analysis complete. {len(df_scores)} predictions generated from {len(successful_scorers)} scorers.")
        
    except Exception as e:
        print(f"Error calling AlphaGenome: {e}")
        results["error"] = str(e)
    
    return results

def process_splice_predictions(df_scores: pd.DataFrame) -> Dict[str, Any]:
    """Process splice-related predictions from AlphaGenome results."""
    splice_types = ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS']
    splice_data = df_scores[df_scores['output_type'].isin(splice_types)]
    
    if splice_data.empty:
        return {"message": "No splice predictions available"}
    
    results = {}
    
    for splice_type in splice_types:
        type_data = splice_data[splice_data['output_type'] == splice_type]
        if not type_data.empty:
            # get significant predictions (high absolute quantile scores)
            significant = type_data[abs(type_data['quantile_score']) > 0.5]
            
            results[splice_type.lower()] = {
                "total_predictions": len(type_data),
                "significant_predictions": len(significant),
                "max_score": type_data['raw_score'].max() if len(type_data) > 0 else None,
                "min_score": type_data['raw_score'].min() if len(type_data) > 0 else None,
                "mean_quantile": type_data['quantile_score'].mean() if len(type_data) > 0 else None,
                "top_tissues": get_top_tissues(type_data) if len(type_data) > 0 else []
            }
    
    return results

def process_expression_predictions(df_scores: pd.DataFrame) -> Dict[str, Any]:
    """Process gene expression predictions from AlphaGenome results."""
    expression_types = ['RNA_SEQ', 'CAGE', 'PROCAP']
    expression_data = df_scores[df_scores['output_type'].isin(expression_types)]
    
    if expression_data.empty:
        return {"message": "No expression predictions available"}
    
    results = {}
    
    for expr_type in expression_types:
        type_data = expression_data[expression_data['output_type'] == expr_type]
        if not type_data.empty:
            # analyze by tissue type
            tissue_analysis = analyze_by_tissue(type_data)
            
            results[expr_type.lower()] = {
                "total_predictions": len(type_data),
                "upregulated_count": len(type_data[type_data['quantile_score'] > 0.2]),
                "downregulated_count": len(type_data[type_data['quantile_score'] < -0.2]),
                "mean_effect": type_data['quantile_score'].mean() if len(type_data) > 0 else None,
                "max_effect": type_data['quantile_score'].max() if len(type_data) > 0 else None,
                "tissue_analysis": tissue_analysis
            }
    
    return results

def process_chromatin_predictions(df_scores: pd.DataFrame) -> Dict[str, Any]:
    """Process chromatin accessibility and histone modification predictions."""
    chromatin_types = ['ATAC', 'DNASE', 'CHIP_HISTONE', 'CHIP_TF']
    chromatin_data = df_scores[df_scores['output_type'].isin(chromatin_types)]
    
    if chromatin_data.empty:
        return {"message": "No chromatin predictions available"}
    
    results = {}
    
    for chrom_type in chromatin_types:
        type_data = chromatin_data[chromatin_data['output_type'] == chrom_type]
        if not type_data.empty:
            results[chrom_type.lower()] = {
                "total_predictions": len(type_data),
                "accessibility_increase": len(type_data[type_data['quantile_score'] > 0.3]),
                "accessibility_decrease": len(type_data[type_data['quantile_score'] < -0.3]),
                "mean_effect": type_data['quantile_score'].mean(),
                "cell_types": get_top_cell_types(type_data)
            }
    
    return results

def get_summary_statistics(df_scores: pd.DataFrame) -> Dict[str, Any]:
    """Get overall summary statistics for the variant predictions."""
    if df_scores.empty:
        return {"total_predictions": 0, "error": "No data available"}
    
    try:
        return {
            "total_predictions": len(df_scores),
            "output_types": df_scores['output_type'].unique().tolist() if 'output_type' in df_scores.columns else [],
            "unique_tissues": df_scores['biosample_name'].nunique() if 'biosample_name' in df_scores.columns else 0,
            "significant_effects": len(df_scores[abs(df_scores['quantile_score']) > 0.5]) if 'quantile_score' in df_scores.columns else 0,
            "mean_quantile_score": df_scores['quantile_score'].mean() if 'quantile_score' in df_scores.columns else 0,
            "max_absolute_effect": abs(df_scores['quantile_score']).max() if 'quantile_score' in df_scores.columns else 0
        }
    except Exception as e:
        return {"error": f"Error computing summary statistics: {e}"}

def get_top_predictions(df_scores: pd.DataFrame, top_n: int = 10) -> List[Dict[str, Any]]:
    """Get the top predictions by absolute quantile score."""
    if df_scores.empty:
        return []
    
    try:
        # get top and bottom scores separately without combining to avoid drop_duplicates issues
        top_scores = df_scores.nlargest(top_n, 'quantile_score', keep='all')
        bottom_scores = df_scores.nsmallest(top_n, 'quantile_score', keep='all')
        
        predictions = []
        
        # add top scores
        for _, row in top_scores.iterrows():
            predictions.append({
                "output_type": row.get('output_type', 'Unknown'),
                "biosample_name": row.get('biosample_name', 'Unknown'),
                "track_name": row.get('track_name', 'Unknown'),
                "raw_score": row.get('raw_score', 0),
                "quantile_score": row.get('quantile_score', 0),
                "effect_direction": "increase" if row.get('quantile_score', 0) > 0 else "decrease",
                "significance": "high" if abs(row.get('quantile_score', 0)) > 0.7 else "moderate" if abs(row.get('quantile_score', 0)) > 0.3 else "low"
            })
        
        # add bottom scores (avoid duplicates manually)
        for _, row in bottom_scores.iterrows():
            duplicate = False
            for existing in predictions:
                if (existing['output_type'] == row.get('output_type', 'Unknown') and 
                    existing['biosample_name'] == row.get('biosample_name', 'Unknown') and
                    existing['quantile_score'] == row.get('quantile_score', 0)):
                    duplicate = True
                    break
            
            if not duplicate:
                predictions.append({
                    "output_type": row.get('output_type', 'Unknown'),
                    "biosample_name": row.get('biosample_name', 'Unknown'),
                    "track_name": row.get('track_name', 'Unknown'),
                    "raw_score": row.get('raw_score', 0),
                    "quantile_score": row.get('quantile_score', 0),
                    "effect_direction": "increase" if row.get('quantile_score', 0) > 0 else "decrease",
                    "significance": "high" if abs(row.get('quantile_score', 0)) > 0.7 else "moderate" if abs(row.get('quantile_score', 0)) > 0.3 else "low"
                })
        
        # sort by absolute quantile score
        predictions.sort(key=lambda x: abs(x['quantile_score']), reverse=True)
        
        return predictions
        
    except Exception as e:
        print(f"Warning: Error in get_top_predictions: {e}")
        return []

def get_top_tissues(data: pd.DataFrame, top_n: int = 5) -> List[Dict[str, Any]]:
    """Get tissues with the strongest effects."""
    if 'biosample_name' not in data.columns:
        return []
    
    tissue_effects = data.groupby('biosample_name')['quantile_score'].agg(['mean', 'max', 'min', 'count']).reset_index()
    tissue_effects['abs_max_effect'] = tissue_effects[['max', 'min']].abs().max(axis=1)
    
    top_tissues = tissue_effects.nlargest(top_n, 'abs_max_effect')
    
    results = []
    for _, row in top_tissues.iterrows():
        results.append({
            "tissue": row['biosample_name'],
            "mean_effect": row['mean'],
            "max_effect": row['max'],
            "min_effect": row['min'],
            "prediction_count": row['count']
        })
    
    return results

def analyze_by_tissue(data: pd.DataFrame) -> Dict[str, Any]:
    """Analyze predictions by tissue/cell type."""
    if 'biosample_name' not in data.columns or data.empty:
        return {}
    
    tissue_stats = data.groupby('biosample_name').agg({
        'quantile_score': ['count', 'mean', 'std', 'max', 'min']
    }).round(3)

    tissue_stats.columns = ['_'.join(col).strip() for col in tissue_stats.columns.values]
    tissue_stats = tissue_stats.reset_index()
    
    return {
        "tissue_count": len(tissue_stats),
        "top_upregulated": tissue_stats.nlargest(3, 'quantile_score_mean')['biosample_name'].tolist(),
        "top_downregulated": tissue_stats.nsmallest(3, 'quantile_score_mean')['biosample_name'].tolist(),
    }

def get_top_cell_types(data: pd.DataFrame, top_n: int = 5) -> List[str]:
    """Get cell types with strongest predicted effects."""
    if 'biosample_name' not in data.columns or data.empty:
        return []
    
    cell_effects = data.groupby('biosample_name')['quantile_score'].apply(lambda x: abs(x).max()).reset_index()
    top_cells = cell_effects.nlargest(top_n, 'quantile_score')
    
    return top_cells['biosample_name'].tolist()

def categorize_alphagenome_effect(quantile_score: float) -> str:
    """Categorize the effect strength based on quantile score."""
    abs_score = abs(quantile_score)
    if abs_score >= 0.8:
        return "very_strong"
    elif abs_score >= 0.5:
        return "strong"
    elif abs_score >= 0.3:
        return "moderate"
    elif abs_score >= 0.1:
        return "weak"
    else:
        return "minimal"

def get_splice_site_analysis(df_scores: pd.DataFrame) -> Dict[str, Any]:
    """Specific analysis for splice site predictions."""
    splice_data = df_scores[df_scores['output_type'] == 'SPLICE_SITES']
    
    if splice_data.empty:
        return {"message": "No splice site predictions available"}
    
    # analyze donor vs acceptor effects if track names indicate this
    donor_effects = splice_data[splice_data['track_name'].str.contains('donor', case=False, na=False)]
    acceptor_effects = splice_data[splice_data['track_name'].str.contains('acceptor', case=False, na=False)]
    
    return {
        "total_splice_predictions": len(splice_data),
        "donor_site_effects": len(donor_effects),
        "acceptor_site_effects": len(acceptor_effects),
        "strongest_effect": splice_data.loc[splice_data['quantile_score'].abs().idxmax()] if len(splice_data) > 0 else None,
        "mean_splice_effect": splice_data['quantile_score'].mean() if len(splice_data) > 0 else None
    }
