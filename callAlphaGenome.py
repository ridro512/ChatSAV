#!/usr/bin/env python3
"""
callAlphaGenome.py - Interface to AlphaGenome API for variant effect prediction

"""

import os
import json
from typing import Dict, Any, Optional, List
import warnings

# Suppress protobuf warnings that are common with AlphaGenome
warnings.filterwarnings("ignore", message="Protobuf gencode version.*is exactly one major version older")


def call_alphagenome(variant_coord: str, hg: str, sequence_length: str = "100KB", 
                    output_types: Optional[List[str]] = None,
                    tissues: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Call AlphaGenome API to predict variant effects.
    
    Args:
        variant_coord: Variant coordinates in format "chr:pos:ref:alt"
        hg: Genome build ("37" or "38") - Note: AlphaGenome primarily supports hg38
        sequence_length: Sequence length for prediction ("2KB", "16KB", "100KB", "500KB", "1MB")
        output_types: List of output types to predict (default: comprehensive set for splice analysis)
        tissues: List of tissue ontology terms (optional)
        
    Returns:
        Dictionary containing AlphaGenome predictions or error message
    """
    try:
        # Import core AlphaGenome modules
        from alphagenome.data import genome
        from alphagenome.models import dna_client
        
        # Try to import variant_scorers
        variant_scorers = None
        try:
            from alphagenome.models import variant_scorers
        except ImportError:
            try:
                from alphagenome.models import variant_scorers
            except ImportError:
                print("Warning: variant_scorers not available. Continuing with basic predictions.")
        
    except ImportError as e:
        return {
            "error": f"AlphaGenome library not installed or incomplete. Please install with: pip install alphagenome. Error: {str(e)}"
        }
    
    # Get API key from environment
    api_key = os.getenv('ALPHA_GENOME_API_KEY')
    if not api_key:
        return {"error": "ALPHA_GENOME_API_KEY environment variable not set"}
    
    # Parse variant coordinates
    try:
        chrom, pos, ref, alt = variant_coord.split(":")
        pos = int(pos)
    except ValueError:
        return {"error": f"Invalid variant format: {variant_coord}. Expected chr:pos:ref:alt"}
    
    # Validate sequence length and convert to AlphaGenome format
    try:
        length_map = {
            "2KB": dna_client.SEQUENCE_LENGTH_2KB,
            "16KB": dna_client.SEQUENCE_LENGTH_16KB,
            "100KB": dna_client.SEQUENCE_LENGTH_100KB,
            "500KB": dna_client.SEQUENCE_LENGTH_500KB,
            "1MB": dna_client.SEQUENCE_LENGTH_1MB
        }
    except AttributeError:
        # Fallback if constants not available
        length_map = {
            "2KB": 2048,
            "16KB": 16384,
            "100KB": 131072,
            "500KB": 524288,
            "1MB": 1048576
        }
    
    if sequence_length not in length_map:
        return {"error": f"Invalid sequence length: {sequence_length}. Must be one of {list(length_map.keys())}"}
    
    sequence_length_bp = length_map[sequence_length]
    
    # Set comprehensive default output types for splice variant analysis
    if output_types is None:
        output_types = ["RNA_SEQ", "ATAC", "DNASE", "CAGE", "PROCAP"]  # Most relevant for splice analysis
    
    # Convert output types to AlphaGenome format
    alphagenome_outputs = []
    try:
        output_type_map = {
            "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
            "ATAC": dna_client.OutputType.ATAC,
            "CAGE": dna_client.OutputType.CAGE,
            "CHIP_TF": dna_client.OutputType.CHIP_TF,
            "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
            "CONTACT_MAPS": dna_client.OutputType.CONTACT_MAPS,
            "DNASE": dna_client.OutputType.DNASE,
            "PROCAP": dna_client.OutputType.PROCAP,
            "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
            "SPLICE_SITE_USAGE": dna_client.OutputType.SPLICE_SITE_USAGE
        }
        
        for output_type in output_types:
            if output_type in output_type_map:
                alphagenome_outputs.append(output_type_map[output_type])
            else:
                return {"error": f"Invalid output type: {output_type}. Must be one of {list(output_type_map.keys())}"}
                
    except AttributeError as e:
        return {"error": f"AlphaGenome OutputType not available: {str(e)}"}
    
    if hg == "37":
        print("Warning: AlphaGenome primarily supports GRCh38/hg38. Results with hg37 coordinates may be less accurate.")
    
    try:
        # Create AlphaGenome model client
        print("Creating AlphaGenome model client...")
        model = dna_client.create(api_key)
        
        # Create variant object
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt
        )
        
        # Create interval centered on variant
        half_length = sequence_length_bp // 2
        start = max(1, pos - half_length)
        end = pos + half_length
        
        interval = genome.Interval(
            chromosome=chrom,
            start=start,
            end=end
        ).resize(sequence_length_bp)
        
        # Set up ontology terms (tissue types)
        ontology_terms = tissues if tissues else []
        
        # Make variant prediction
        print(f"Making AlphaGenome prediction for {variant_coord}...")
        print(f"Interval: {interval}")
        print(f"Sequence length: {sequence_length} ({sequence_length_bp:,} bp)")
        print(f"Output types: {', '.join(output_types)}")
        
        variant_output = model.predict_variant(
            interval=interval,
            variant=variant,
            requested_outputs=alphagenome_outputs,
            ontology_terms=ontology_terms
        )
        
        # Get variant scores using recommended and splice-specific scorers
        variant_scores_results = {}
        if variant_scorers is not None:
            try:
                print("Computing variant scores with multiple scorers...")
                
                # Get recommended scorers for each output type
                for output_type in output_types:
                    scorers_to_use = []
                    
                    # Try to get recommended scorer
                    try:
                        recommended_scorers = variant_scorers.get_recommended_scorers('human')
                        for scorer in recommended_scorers:
                            if hasattr(scorer, 'requested_output') and str(scorer.requested_output).endswith(output_type):
                                scorers_to_use.append(scorer)
                    except Exception as e:
                        print(f"Warning: Could not get recommended scorers: {e}")
                    
                    # Add splice-specific scorers for RNA_SEQ
                    if output_type == "RNA_SEQ":
                        try:
                            # Add splice junction scorer
                            splice_scorer = variant_scorers.SpliceJunctionScorer()
                            scorers_to_use.append(splice_scorer)
                            
                            
                        except Exception as e:
                            print(f"Warning: Could not add splice-specific scorers: {e}")
                    
                    # Add splice-specific scorers for SPLICE_SITES and SPLICE_SITE_USAGE
                    if output_type in ["SPLICE_SITES", "SPLICE_SITE_USAGE"]:
                        try:
                            # Add gene mask splicing scorer (only for splice-specific outputs)
                            gene_splice_scorer = variant_scorers.GeneMaskSplicingScorer(
                                requested_output=output_type_map[output_type],
                                width=10001
                            )
                            scorers_to_use.append(gene_splice_scorer)
                            
                        except Exception as e:
                            print(f"Warning: Could not add gene mask splicing scorer: {e}")
                    
                    # Add center mask scorer for comprehensive analysis
                    try:
                        # Import AggregationType
                        aggregation_type = variant_scorers.AggregationType.DIFF_MEAN  # Default to mean aggregation
                        
                        center_scorer = variant_scorers.CenterMaskScorer(
                            requested_output=output_type_map[output_type],
                            width=10001,
                            aggregation_type=aggregation_type
                        )
                        scorers_to_use.append(center_scorer)
                    except Exception as e:
                        print(f"Warning: Could not add center mask scorer: {e}")
                    
                    # Compute scores if we have scorers
                    if scorers_to_use:
                        try:
                            scores = model.score_variant(
                                interval=interval,
                                variant=variant,
                                variant_scorers=scorers_to_use
                            )
                            if scores:
                                variant_scores_results[output_type] = scores
                                print(f"✓ Computed {len(scores)} variant scores for {output_type}")
                        except Exception as e:
                            print(f"Warning: Could not compute variant scores for {output_type}: {e}")
                
            except Exception as e:
                print(f"Warning: Variant scoring failed: {e}")
        else:
            print("Note: Variant scoring not available - using prediction data only")
        
        # Process and format results
        result = {
            "variant": variant_coord,
            "genome_build": f"hg{hg}",
            "sequence_length": sequence_length,
            "interval": {
                "chromosome": interval.chromosome,
                "start": interval.start,
                "end": interval.end,
                "length": interval.end - interval.start
            },
            "predictions": {},
            "variant_scores": {}
        }
        
        # Process predictions for each output type
        for output_type in output_types:
            output_key = output_type.lower()
            
            prediction_data = {
                "reference_available": False,
                "alternate_available": False,
                "reference_summary": {},
                "alternate_summary": {},
                "prediction_diff": {}
            }
            
            # Get reference and alternate predictions
            try:
                if hasattr(variant_output, 'reference'):
                    ref_output = getattr(variant_output.reference, output_key, None)
                    if ref_output is not None:
                        prediction_data["reference_available"] = True
                        prediction_data["reference_summary"] = summarize_track_data(ref_output)
                
                if hasattr(variant_output, 'alternate'):
                    alt_output = getattr(variant_output.alternate, output_key, None)
                    if alt_output is not None:
                        prediction_data["alternate_available"] = True
                        prediction_data["alternate_summary"] = summarize_track_data(alt_output)
                
                # Calculate basic difference if both available
                if prediction_data["reference_available"] and prediction_data["alternate_available"]:
                    ref_stats = prediction_data["reference_summary"].get("value_stats", {})
                    alt_stats = prediction_data["alternate_summary"].get("value_stats", {})
                    
                    if ref_stats and alt_stats:
                        prediction_data["prediction_diff"] = {
                            "mean_diff": alt_stats.get("mean", 0) - ref_stats.get("mean", 0),
                            "max_diff": alt_stats.get("max", 0) - ref_stats.get("max", 0),
                            "effect_magnitude": abs(alt_stats.get("mean", 0) - ref_stats.get("mean", 0))
                        }
                        
            except Exception as e:
                prediction_data["error"] = f"Error processing {output_type}: {str(e)}"
            
            result["predictions"][output_type] = prediction_data
        
        # Process variant scores
        for output_type, scores_list in variant_scores_results.items():
            result["variant_scores"][output_type] = []
            
            for i, scores_data in enumerate(scores_list):
                processed_scores = process_variant_scores(scores_data, output_type, i)
                result["variant_scores"][output_type].append(processed_scores)
        
        return result
        
    except Exception as e:
        return {"error": f"AlphaGenome prediction failed: {str(e)}"}


def summarize_track_data(track_data) -> Dict[str, Any]:
    """
    Summarize track data from AlphaGenome predictions.
    
    Args:
        track_data: Track data from AlphaGenome prediction
        
    Returns:
        Dictionary with summary statistics
    """
    summary = {
        "tracks_available": False,
        "track_count": 0,
        "value_stats": {},
        "data_shape": None
    }
    
    try:
        if hasattr(track_data, 'values') and track_data.values is not None:
            values = track_data.values
            summary["tracks_available"] = True
            summary["data_shape"] = list(values.shape)
            summary["track_count"] = values.shape[1] if len(values.shape) > 1 else 1
            
            # Calculate basic statistics
            try:
                import numpy as np
                
                # Flatten for overall stats - with better empty array handling
                if hasattr(values, 'flatten'):
                    flat_values = values.flatten()
                else:
                    flat_values = values
                
                if len(flat_values) == 0:
                    summary["value_stats"] = {"error": "No values available"}
                else:
                    # Filter out NaN and infinite values before calculations
                    valid_values = flat_values[np.isfinite(flat_values)]
                    
                    if len(valid_values) == 0:
                        summary["value_stats"] = {"error": "No valid (finite) values available"}
                    else:
                        summary["value_stats"] = {
                            "mean": float(np.mean(valid_values)),
                            "std": float(np.std(valid_values)),
                            "min": float(np.min(valid_values)),
                            "max": float(np.max(valid_values)),
                            "non_zero_count": int(np.count_nonzero(valid_values)),
                            "total_count": int(len(valid_values)),
                            "valid_ratio": len(valid_values) / len(flat_values)
                        }
                        
                        # Add percentiles only if we have enough valid values
                        if len(valid_values) >= 4:
                            summary["value_stats"]["percentiles"] = {
                                "25th": float(np.percentile(valid_values, 25)),
                                "50th": float(np.percentile(valid_values, 50)),
                                "75th": float(np.percentile(valid_values, 75)),
                                "95th": float(np.percentile(valid_values, 95))
                            }
                
            except ImportError:
                summary["value_stats"] = {"error": "numpy not available for statistics"}
            except Exception as e:
                summary["value_stats"] = {"error": f"Statistics calculation failed: {str(e)}"}
                
        if hasattr(track_data, 'interval'):
            summary["interval"] = {
                "chromosome": track_data.interval.chromosome,
                "start": track_data.interval.start,
                "end": track_data.interval.end,
                "length": track_data.interval.end - track_data.interval.start
            }
            
        # Try to get metadata if available
        if hasattr(track_data, 'metadata'):
            summary["metadata"] = track_data.metadata
            
    except Exception as e:
        summary["error"] = f"Failed to summarize track data: {str(e)}"
    
    return summary


def process_variant_scores(scores_data, output_type: str, scorer_index: int) -> Dict[str, Any]:
    """
    Process variant scores from AlphaGenome into a readable format.
    
    Args:
        scores_data: AnnData object containing variant scores
        output_type: The output type being processed
        scorer_index: Index of the scorer used
        
    Returns:
        Dictionary with processed score information
    """
    try:
        # Import required libraries
        try:
            import pandas as pd
            import numpy as np
        except ImportError:
            return {"error": "pandas and numpy required for score processing"}
        
        result = {
            "output_type": output_type,
            "scorer_index": scorer_index,
            "scorer_type": "unknown",
            "total_genes": 0,
            "total_tracks": 0,
            "top_affected_genes": [],
            "summary_stats": {},
            "splice_specific_results": {}
        }
        
        # Get scorer information
        if hasattr(scores_data, 'uns') and 'variant_scorer' in scores_data.uns:
            scorer_info = scores_data.uns['variant_scorer']
            result["scorer_type"] = type(scorer_info).__name__
        
        if hasattr(scores_data, 'X') and scores_data.X is not None:
            # Basic shape information
            result["total_genes"] = scores_data.X.shape[0]
            result["total_tracks"] = scores_data.X.shape[1]
            
            # Gene information
            if hasattr(scores_data, 'obs') and not scores_data.obs.empty:
                genes_df = scores_data.obs.copy()
                
                # Calculate mean absolute effect score per gene with error handling
                if scores_data.X.size > 0:
                    mean_scores = np.abs(scores_data.X).mean(axis=1)
                    max_scores = np.abs(scores_data.X).max(axis=1)
                    
                    # Handle NaN values
                    mean_scores = np.nan_to_num(mean_scores, nan=0.0)
                    max_scores = np.nan_to_num(max_scores, nan=0.0)
                    
                    genes_df['mean_abs_score'] = mean_scores
                    genes_df['max_abs_score'] = max_scores
                    
                    # Sort by maximum absolute score and get top genes (only if we have valid scores)
                    if np.any(max_scores > 0):
                        top_genes = genes_df.nlargest(10, 'max_abs_score')
                        
                        result["top_affected_genes"] = []
                        for idx, row in top_genes.iterrows():
                            gene_info = {
                                "gene_id": row.get('gene_id', 'N/A'),
                                "gene_name": row.get('gene_name', 'N/A'),
                                "gene_type": row.get('gene_type', 'N/A'),
                                "strand": row.get('strand', 'N/A'),
                                "mean_abs_score": float(row['mean_abs_score']),
                                "max_abs_score": float(row['max_abs_score'])
                            }
                            result["top_affected_genes"].append(gene_info)
            
            # Overall summary statistics with better error handling
            if scores_data.X.size > 0:
                all_scores = scores_data.X.flatten()
                # Remove NaN and infinite values
                valid_scores = all_scores[np.isfinite(all_scores)]
                
                if len(valid_scores) > 0:
                    result["summary_stats"] = {
                        "mean_score": float(np.mean(valid_scores)),
                        "std_score": float(np.std(valid_scores)),
                        "min_score": float(np.min(valid_scores)),
                        "max_score": float(np.max(valid_scores)),
                        "significant_scores": int(np.sum(np.abs(valid_scores) > 0.1)),
                        "total_scores": len(valid_scores),
                        "valid_score_ratio": len(valid_scores) / len(all_scores)
                    }
                else:
                    result["summary_stats"] = {
                        "mean_score": 0.0,
                        "std_score": 0.0,
                        "min_score": 0.0,
                        "max_score": 0.0,
                        "significant_scores": 0,
                        "total_scores": 0,
                        "valid_score_ratio": 0.0,
                        "note": "No valid scores found"
                    }
            else:
                result["summary_stats"] = {
                    "error": "No score data available"
                }
            
            # Add splice-specific analysis for RNA_SEQ with splice scorers
            if output_type == "RNA_SEQ" and result["scorer_type"] in ["SpliceJunctionScorer", "GeneMaskSplicingScorer"]:
                result["splice_specific_results"] = analyze_splice_scores(scores_data)
        
        return result
        
    except Exception as e:
        return {"error": f"Failed to process variant scores: {str(e)}"}


def analyze_splice_scores(scores_data) -> Dict[str, Any]:
    """
    Analyze splice-specific variant scores.
    
    Args:
        scores_data: AnnData object with splice variant scores
        
    Returns:
        Dictionary with splice-specific analysis
    """
    try:
        import numpy as np
        
        splice_analysis = {
            "splice_disruption_genes": [],
            "splice_enhancement_genes": [],
            "total_splice_effects": 0,
            "strong_splice_effects": 0
        }
        
        if hasattr(scores_data, 'X') and scores_data.X is not None and scores_data.X.size > 0:
            # Threshold for significant splice effects
            splice_threshold = 0.2
            strong_threshold = 0.5
            
            # Analyze per gene
            if hasattr(scores_data, 'obs') and not scores_data.obs.empty:
                for idx, row in scores_data.obs.reset_index(drop=True).iterrows():
                    if not isinstance(idx, int) or idx >= scores_data.X.shape[0]:
                        continue
                        
                    gene_scores = scores_data.X[idx, :]
                    
                    # Handle empty or invalid scores
                    if gene_scores.size == 0:
                        continue
                    
                    # Remove NaN and infinite values
                    valid_scores = gene_scores[np.isfinite(gene_scores)]
                    if len(valid_scores) == 0:
                        continue
                    
                    max_score = np.max(valid_scores)
                    min_score = np.min(valid_scores)
                    
                    # Check for splice disruption (negative scores typically indicate disruption)
                    if min_score < -splice_threshold:
                        splice_analysis["splice_disruption_genes"].append({
                            "gene_name": row.get('gene_name', 'Unknown'),
                            "gene_id": row.get('gene_id', 'Unknown'),
                            "disruption_score": float(min_score),
                            "strength": "strong" if min_score < -strong_threshold else "moderate"
                        })
                    
                    # Check for splice enhancement (positive scores)
                    if max_score > splice_threshold:
                        splice_analysis["splice_enhancement_genes"].append({
                            "gene_name": row.get('gene_name', 'Unknown'),
                            "gene_id": row.get('gene_id', 'Unknown'),
                            "enhancement_score": float(max_score),
                            "strength": "strong" if max_score > strong_threshold else "moderate"
                        })
                
                # Count effects across all valid scores
                all_scores = scores_data.X.flatten()
                valid_all_scores = all_scores[np.isfinite(all_scores)]
                
                if len(valid_all_scores) > 0:
                    splice_analysis["total_splice_effects"] = int(np.sum(np.abs(valid_all_scores) > splice_threshold))
                    splice_analysis["strong_splice_effects"] = int(np.sum(np.abs(valid_all_scores) > strong_threshold))
        
        return splice_analysis
        
    except Exception as e:
        return {"error": f"Splice analysis failed: {str(e)}"}


def get_available_output_types() -> List[str]:
    """
    Get list of available AlphaGenome output types.
    
    Returns:
        List of available output type names
    """
    return ["RNA_SEQ", "ATAC", "CAGE", "CHIP_TF", "CHIP_HISTONE", "CONTACT_MAPS", "DNASE", "PROCAP", "SPLICE_SITES", "SPLICE_SITE_USAGE"]


def get_sequence_length_info() -> Dict[str, Dict[str, Any]]:
    """
    Get information about available sequence lengths.
    
    Returns:
        Dictionary with sequence length information
    """
    try:
        from alphagenome.models import dna_client
        
        return {
            "2KB": {
                "base_pairs": getattr(dna_client, 'SEQUENCE_LENGTH_2KB', 2048),
                "description": "Short range - local splice sites and immediate regulatory elements",
                "use_case": "Point mutations with local effects"
            },
            "16KB": {
                "base_pairs": getattr(dna_client, 'SEQUENCE_LENGTH_16KB', 16384),
                "description": "Medium range - nearby enhancers and regulatory regions",
                "use_case": "Variants affecting local gene regulation"
            },
            "100KB": {
                "base_pairs": getattr(dna_client, 'SEQUENCE_LENGTH_100KB', 131072),
                "description": "Standard range - gene and proximal regulatory elements",
                "use_case": "Most variants - balanced resolution and context"
            },
            "500KB": {
                "base_pairs": getattr(dna_client, 'SEQUENCE_LENGTH_500KB', 524288),
                "description": "Long range - distal regulatory elements and TAD boundaries",
                "use_case": "Variants in regulatory deserts or structural variants"
            },
            "1MB": {
                "base_pairs": getattr(dna_client, 'SEQUENCE_LENGTH_1MB', 1048576),
                "description": "Very long range - chromosome topology and long-range interactions",
                "use_case": "Large structural variants or complex regulatory landscapes"
            }
        }
    except ImportError:
        return {
            "2KB": {"base_pairs": 2048, "description": "Short range analysis"},
            "16KB": {"base_pairs": 16384, "description": "Medium range analysis"},
            "100KB": {"base_pairs": 131072, "description": "Standard range analysis"},
            "500KB": {"base_pairs": 524288, "description": "Long range analysis"},
            "1MB": {"base_pairs": 1048576, "description": "Very long range analysis"}
        }


def get_common_tissue_ontologies() -> Dict[str, str]:
    """
    Get mapping of common tissue names to ontology terms.
    
    Returns:
        Dictionary mapping tissue names to UBERON ontology terms
    """
    return {
        "brain": "UBERON:0000955",
        "heart": "UBERON:0000948",
        "liver": "UBERON:0002107",
        "lung": "UBERON:0002048",
        "kidney": "UBERON:0002113",
        "muscle": "UBERON:0000993",
        "colon": "UBERON:0001157",
        "blood": "UBERON:0000178",
        "skin": "UBERON:0002097",
        "pancreas": "UBERON:0001264",
        "stomach": "UBERON:0000945",
        "breast": "UBERON:0000310",
        "prostate": "UBERON:0002367",
        "esophagus": "UBERON:0001043"
    }
