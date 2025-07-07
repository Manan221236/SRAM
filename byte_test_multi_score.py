#!/usr/bin/env python3
"""
Multi-Metric Cryptographic Attack Analysis
Analyzes results from 10-iteration multi-metric attack testbench
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cosine as cosine_distance
from sklearn.metrics import mean_squared_error

warnings.filterwarnings("ignore")

# Configuration
W, H = 128, 128
TOTAL = W * H
NUM_ITERATIONS = 10

# Color palette for metrics
METRIC_COLORS = {
    'pearson': '#1f77b4',    # Blue
    'spearman': '#ff7f0e',   # Orange
    'cosine': '#2ca02c',     # Green
    'mse': '#d62728',        # Red
    'composite': '#9467bd'   # Purple
}

def load_img(path: Path) -> np.ndarray:
    """Load image from hex file, handling missing/corrupted data"""
    pix = []
    try:
        with path.open() as fh:
            for ln in fh:
                t = ln.strip()
                if not t or t.startswith("#"):
                    continue
                if t == "XXXXXXXXXXXXXXX":
                    pix.append(0)
                else:
                    try:
                        pix.append(int(t, 16) & 0xFF)
                    except ValueError:
                        pix.append(0)
    except FileNotFoundError:
        print(f"Warning: {path} not found, returning empty array")
        return np.array([])
    
    return np.array(pix[:TOTAL], dtype=np.uint8)

def calculate_all_metrics(img, ref):
    """Calculate all 5 metrics: Pearson, Spearman, Cosine, MSE, Composite"""
    if len(img) == 0 or len(ref) == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    
    # Ensure same length
    img = img[:TOTAL]
    ref = ref[:TOTAL]
    
    # Pearson correlation
    try:
        pearson_corr, _ = pearsonr(img, ref)
        if np.isnan(pearson_corr):
            pearson_corr = 0.0
    except:
        pearson_corr = 0.0
    
    # Spearman rank correlation
    try:
        spearman_corr, _ = spearmanr(img, ref)
        if np.isnan(spearman_corr):
            spearman_corr = 0.0
    except:
        spearman_corr = 0.0
    
    # Cosine similarity
    try:
        # cosine_distance returns distance, so similarity = 1 - distance
        cosine_sim = 1.0 - cosine_distance(img, ref)
        if np.isnan(cosine_sim):
            cosine_sim = 0.0
    except:
        cosine_sim = 0.0
    
    # MSE (inverted to similarity)
    try:
        mse = mean_squared_error(ref, img)
        mse_sim = 1.0 / (1.0 + (mse / 65025.0))  # Normalize by 255^2
    except:
        mse_sim = 0.0
    
    # Composite score (using same weights as Verilog)
    # Normalize correlations to [0,1]
    pearson_norm = (pearson_corr + 1.0) / 2.0
    spearman_norm = (spearman_corr + 1.0) / 2.0
    cosine_norm = (cosine_sim + 1.0) / 2.0
    
    composite = 0.40 * pearson_norm + 0.25 * spearman_norm + 0.20 * cosine_norm + 0.15 * mse_sim
    
    return pearson_corr, spearman_corr, cosine_sim, mse_sim, composite

def load_summary_data():
    """Load summary data from Verilog-generated file"""
    summary_file = Path("multi_metric_summary.txt")
    if not summary_file.exists():
        print("Warning: multi_metric_summary.txt not found")
        return None
    
    data = {
        'keys': [],
        'pearson': [],
        'spearman': [],
        'cosine': [],
        'mse': [],
        'composite': []
    }
    
    try:
        with summary_file.open() as f:
            lines = f.readlines()
            for line in lines:
                if line.strip() and '|' in line and not line.startswith('Iter'):
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 7:
                        data['keys'].append(parts[1])
                        data['pearson'].append(float(parts[2]))
                        data['spearman'].append(float(parts[3]))
                        data['cosine'].append(float(parts[4]))
                        data['mse'].append(float(parts[5]))
                        data['composite'].append(float(parts[6]))
    except Exception as e:
        print(f"Error parsing summary file: {e}")
        return None
    
    return data

def create_metric_comparison_plots(iteration_metrics, ideal_metrics):
    """Create comprehensive metric comparison visualizations"""
    fig = plt.figure(figsize=(20, 12))
    
    # Plot 1: All metrics over iterations
    plt.subplot(2, 4, 1)
    iterations = range(NUM_ITERATIONS)
    
    for metric in ['pearson', 'spearman', 'cosine', 'mse', 'composite']:
        values = iteration_metrics[metric]
        plt.plot(iterations, values, 'o-', linewidth=2, markersize=6, 
                color=METRIC_COLORS[metric], label=metric.title())
        
        # Add ideal line if available
        if ideal_metrics and metric in ideal_metrics:
            plt.axhline(y=ideal_metrics[metric], color=METRIC_COLORS[metric], 
                       linestyle='--', alpha=0.5)
    
    plt.xlabel('Iteration')
    plt.ylabel('Metric Value')
    plt.title('All Metrics Over Iterations')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.xticks(iterations)
    
    # Plot 2: Composite score focus
    plt.subplot(2, 4, 2)
    composite_values = iteration_metrics['composite']
    bars = plt.bar(iterations, composite_values, alpha=0.7, 
                   color=METRIC_COLORS['composite'], edgecolor='navy')
    
    if ideal_metrics and 'composite' in ideal_metrics:
        plt.axhline(y=ideal_metrics['composite'], color='red', linestyle='--', 
                   linewidth=2, label=f'Ideal ({ideal_metrics["composite"]:.4f})')
        plt.legend()
    
    plt.xlabel('Iteration')
    plt.ylabel('Composite Score')
    plt.title('Composite Score Progress')
    plt.grid(True, alpha=0.3)
    plt.xticks(iterations)
    
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, composite_values)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{val:.4f}', ha='center', fontsize=9, rotation=90)
    
    # Plot 3: Metric correlation heatmap
    plt.subplot(2, 4, 3)
    metric_data = np.array([iteration_metrics[m] for m in ['pearson', 'spearman', 'cosine', 'mse', 'composite']])
    correlation_matrix = np.corrcoef(metric_data)
    
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0,
                xticklabels=['Pearson', 'Spearman', 'Cosine', 'MSE', 'Composite'],
                yticklabels=['Pearson', 'Spearman', 'Cosine', 'MSE', 'Composite'],
                cbar_kws={'label': 'Correlation'})
    plt.title('Inter-Metric Correlations')
    
    # Plot 4: Improvement analysis
    plt.subplot(2, 4, 4)
    if len(composite_values) > 1:
        improvements = [composite_values[i+1] - composite_values[i] 
                       for i in range(len(composite_values)-1)]
        colors = ['green' if imp > 0 else 'red' for imp in improvements]
        plt.bar(range(1, NUM_ITERATIONS), improvements, alpha=0.7, color=colors)
        plt.xlabel('Iteration Transition')
        plt.ylabel('Composite Score Change')
        plt.title('Iteration-to-Iteration Improvement')
        plt.grid(True, alpha=0.3)
        plt.xticks(range(1, NUM_ITERATIONS), [f'{i-1}→{i}' for i in range(1, NUM_ITERATIONS)])
    
    # Plot 5: Individual metric trends
    plt.subplot(2, 4, 5)
    for metric in ['pearson', 'spearman', 'cosine', 'mse']:
        values = iteration_metrics[metric]
        plt.plot(iterations, values, 'o-', linewidth=2, markersize=4, 
                color=METRIC_COLORS[metric], label=metric.title(), alpha=0.8)
    
    plt.xlabel('Iteration')
    plt.ylabel('Metric Value')
    plt.title('Individual Metric Trends')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xticks(iterations)
    
    # Plot 6: Best vs Worst iteration comparison
    plt.subplot(2, 4, 6)
    best_idx = np.argmax(composite_values)
    worst_idx = np.argmin(composite_values)
    
    metrics = ['Pearson', 'Spearman', 'Cosine', 'MSE', 'Composite']
    best_values = [iteration_metrics['pearson'][best_idx], iteration_metrics['spearman'][best_idx],
                  iteration_metrics['cosine'][best_idx], iteration_metrics['mse'][best_idx],
                  iteration_metrics['composite'][best_idx]]
    worst_values = [iteration_metrics['pearson'][worst_idx], iteration_metrics['spearman'][worst_idx],
                   iteration_metrics['cosine'][worst_idx], iteration_metrics['mse'][worst_idx],
                   iteration_metrics['composite'][worst_idx]]
    
    x = np.arange(len(metrics))
    width = 0.35
    
    plt.bar(x - width/2, best_values, width, label=f'Best (Iter {best_idx})', 
           color='lightgreen', alpha=0.8)
    plt.bar(x + width/2, worst_values, width, label=f'Worst (Iter {worst_idx})', 
           color='lightcoral', alpha=0.8)
    
    plt.xlabel('Metrics')
    plt.ylabel('Values')
    plt.title('Best vs Worst Iteration')
    plt.xticks(x, metrics)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 7: Attack effectiveness
    plt.subplot(2, 4, 7)
    if ideal_metrics:
        effectiveness = []
        for metric in ['pearson', 'spearman', 'cosine', 'mse', 'composite']:
            if metric in ideal_metrics and ideal_metrics[metric] > 0:
                best_val = max(iteration_metrics[metric])
                eff = (best_val / ideal_metrics[metric]) * 100
                effectiveness.append(eff)
            else:
                effectiveness.append(0)
        
        bars = plt.bar(['Pearson', 'Spearman', 'Cosine', 'MSE', 'Composite'], 
                      effectiveness, color=[METRIC_COLORS[m] for m in ['pearson', 'spearman', 'cosine', 'mse', 'composite']], 
                      alpha=0.7)
        plt.xlabel('Metrics')
        plt.ylabel('Effectiveness (%)')
        plt.title('Attack Effectiveness vs Ideal')
        plt.grid(True, alpha=0.3)
        
        # Add percentage labels
        for bar, eff in zip(bars, effectiveness):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                    f'{eff:.1f}%', ha='center', fontsize=10)
    
    # Plot 8: Convergence analysis
    plt.subplot(2, 4, 8)
    plt.plot(iterations, composite_values, 'o-', linewidth=3, markersize=8, 
            color=METRIC_COLORS['composite'], label='Composite Score')
    
    # Fit trend line
    if len(composite_values) > 1:
        z = np.polyfit(iterations, composite_values, 1)
        p = np.poly1d(z)
        plt.plot(iterations, p(iterations), "--", alpha=0.8, color='red', 
                label=f'Trend (slope: {z[0]:+.4f})')
    
    plt.xlabel('Iteration')
    plt.ylabel('Composite Score')
    plt.title('Convergence Analysis')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xticks(iterations)
    
    plt.tight_layout()
    plt.savefig("multi_metric_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()

def create_image_gallery(iteration_images, ref_image, ideal_image=None):
    """Create visual gallery of all iteration results"""
    n_cols = 4
    n_rows = 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 12))
    axes = axes.ravel()
    
    # Show iteration results
    for i in range(min(NUM_ITERATIONS, 10)):
        if i < len(iteration_images):
            axes[i].imshow(iteration_images[i].reshape(H, W), cmap="gray", vmin=0, vmax=255)
            axes[i].set_title(f"Iteration {i}", fontsize=10)
        else:
            axes[i].text(0.5, 0.5, "No Data", ha='center', va='center', 
                        transform=axes[i].transAxes)
        axes[i].axis("off")
    
    # Show reference image
    if len(axes) > 10:
        axes[10].imshow(ref_image.reshape(H, W), cmap="gray", vmin=0, vmax=255)
        axes[10].set_title("Reference\n(Original)", fontsize=10, color='black', weight='bold')
        axes[10].axis("off")
    
    # Show ideal case
    if len(axes) > 11 and ideal_image is not None:
        axes[11].imshow(ideal_image.reshape(H, W), cmap="gray", vmin=0, vmax=255)
        axes[11].set_title("Ideal Case", fontsize=10, color='red', weight='bold')
        axes[11].axis("off")
    elif len(axes) > 11:
        axes[11].text(0.5, 0.5, "Ideal Case\nNot Available", ha='center', va='center', 
                     transform=axes[11].transAxes)
        axes[11].axis("off")
    
    fig.suptitle("Multi-Metric Attack Results Gallery", weight="bold", fontsize=14)
    fig.tight_layout()
    plt.savefig("multi_metric_gallery.png", dpi=300, bbox_inches="tight")
    plt.close()

def generate_detailed_report(iteration_metrics, ideal_metrics, summary_data):
    """Generate comprehensive text report"""
    with open("multi_metric_detailed_report.txt", "w") as f:
        f.write("Multi-Metric Cryptographic Attack Analysis Report\n")
        f.write("=" * 55 + "\n\n")
        
        # Executive Summary
        f.write("EXECUTIVE SUMMARY:\n")
        f.write("-" * 18 + "\n")
        f.write(f"Total iterations analyzed: {NUM_ITERATIONS}\n")
        f.write(f"Image dimensions: {W}×{H} = {TOTAL} pixels\n")
        f.write(f"Scoring method: Multi-metric composite (Pearson 40%, Spearman 25%, Cosine 20%, MSE 15%)\n\n")
        
        # Best Results
        best_idx = np.argmax(iteration_metrics['composite'])
        f.write("BEST ITERATION RESULTS:\n")
        f.write("-" * 24 + "\n")
        f.write(f"Best iteration: {best_idx}\n")
        f.write(f"Composite score: {iteration_metrics['composite'][best_idx]:.6f}\n")
        f.write(f"Pearson correlation: {iteration_metrics['pearson'][best_idx]:.6f}\n")
        f.write(f"Spearman correlation: {iteration_metrics['spearman'][best_idx]:.6f}\n")
        f.write(f"Cosine similarity: {iteration_metrics['cosine'][best_idx]:.6f}\n")
        f.write(f"MSE similarity: {iteration_metrics['mse'][best_idx]:.6f}\n\n")
        
        # Statistical Analysis
        f.write("STATISTICAL ANALYSIS:\n")
        f.write("-" * 21 + "\n")
        for metric in ['composite', 'pearson', 'spearman', 'cosine', 'mse']:
            values = iteration_metrics[metric]
            f.write(f"{metric.upper()} STATISTICS:\n")
            f.write(f"  Mean: {np.mean(values):.6f}\n")
            f.write(f"  Std:  {np.std(values):.6f}\n")
            f.write(f"  Min:  {np.min(values):.6f} (Iteration {np.argmin(values)})\n")
            f.write(f"  Max:  {np.max(values):.6f} (Iteration {np.argmax(values)})\n\n")
        
        # Convergence Analysis
        if len(iteration_metrics['composite']) > 1:
            f.write("CONVERGENCE ANALYSIS:\n")
            f.write("-" * 21 + "\n")
            
            # Overall trend
            composite_values = iteration_metrics['composite']
            trend_slope = np.polyfit(range(len(composite_values)), composite_values, 1)[0]
            
            if trend_slope > 0.001:
                convergence = "IMPROVING"
            elif trend_slope < -0.001:
                convergence = "DEGRADING"
            else:
                convergence = "STABLE"
            
            f.write(f"Overall trend: {convergence} (slope: {trend_slope:+.6f})\n")
            
            # Total improvement
            total_improvement = composite_values[-1] - composite_values[0]
            f.write(f"Total improvement: {total_improvement:+.6f}\n")
            
            # Best improvement step
            improvements = [composite_values[i+1] - composite_values[i] 
                           for i in range(len(composite_values)-1)]
            best_step = np.argmax(improvements)
            f.write(f"Best improvement: Iteration {best_step} → {best_step+1} ({improvements[best_step]:+.6f})\n\n")
        
        # Effectiveness vs Ideal
        if ideal_metrics:
            f.write("EFFECTIVENESS vs IDEAL CASE:\n")
            f.write("-" * 29 + "\n")
            for metric in ['pearson', 'spearman', 'cosine', 'mse', 'composite']:
                if metric in ideal_metrics and ideal_metrics[metric] > 0:
                    best_val = max(iteration_metrics[metric])
                    effectiveness = (best_val / ideal_metrics[metric]) * 100
                    f.write(f"{metric.title():10s}: {effectiveness:6.2f}% of theoretical maximum\n")
            f.write("\n")
        
        # Inter-metric Correlations
        f.write("INTER-METRIC CORRELATIONS:\n")
        f.write("-" * 26 + "\n")
        metrics = ['pearson', 'spearman', 'cosine', 'mse', 'composite']
        metric_data = np.array([iteration_metrics[m] for m in metrics])
        corr_matrix = np.corrcoef(metric_data)
        
        f.write("         " + "".join([f"{m:>9s}" for m in metrics]) + "\n")
        for i, metric in enumerate(metrics):
            f.write(f"{metric:>8s} ")
            for j in range(len(metrics)):
                f.write(f"{corr_matrix[i,j]:8.3f} ")
            f.write("\n")
        f.write("\n")
        
        # Recommendations
        f.write("RECOMMENDATIONS:\n")
        f.write("-" * 16 + "\n")
        
        best_composite = max(iteration_metrics['composite'])
        if best_composite < 0.3:
            f.write("• Attack showing poor results - consider different approach or longer iterations\n")
        elif best_composite < 0.6:
            f.write("• Moderate success - algorithm working but may need parameter tuning\n")
        elif best_composite < 0.8:
            f.write("• Good performance - attack is effective with room for optimization\n")
        else:
            f.write("• Excellent performance - attack is highly effective\n")
        
        if len(iteration_metrics['composite']) > 1:
            if iteration_metrics['composite'][-1] > iteration_metrics['composite'][0]:
                f.write("• Iterations are helping - continue refinement approach\n")
            else:
                f.write("• No improvement from iterations - check algorithm parameters\n")
        
        # Check metric agreement
        metric_correlations = []
        for i in range(len(metrics)-1):
            for j in range(i+1, len(metrics)-1):  # Exclude composite
                metric_correlations.append(corr_matrix[i,j])
        
        avg_corr = np.mean(metric_correlations)
        if avg_corr > 0.8:
            f.write("• High metric agreement - all metrics pointing in same direction\n")
        elif avg_corr < 0.3:
            f.write("• Low metric agreement - consider adjusting composite weights\n")
        
        f.write("\nEnd of Report\n")

def main():
    """Main analysis function"""
    print("Starting Multi-Metric Attack Analysis...")
    print("=" * 50)
    
    # Load reference image
    print("Loading reference images...")
    ref_image = load_img(Path("reference_correct_decode.txt"))
    if ref_image.size == 0:
        print("ERROR: reference_correct_decode.txt missing or empty!")
        return
    
    # Load ideal case
    ideal_image = load_img(Path("correct_trnga_wrong_trngd.txt"))
    
    # Load iteration results
    print("Loading iteration results...")
    iteration_images = []
    iteration_metrics = {
        'pearson': [],
        'spearman': [],
        'cosine': [],
        'mse': [],
        'composite': []
    }
    
    for iter_idx in range(NUM_ITERATIONS):
        img_path = Path(f"final_iteration_{iter_idx}.txt")
        img = load_img(img_path)
        
        if img.size > 0:
            pearson, spearman, cosine, mse, composite = calculate_all_metrics(img, ref_image)
            iteration_images.append(img)
            iteration_metrics['pearson'].append(pearson)
            iteration_metrics['spearman'].append(spearman)
            iteration_metrics['cosine'].append(cosine)
            iteration_metrics['mse'].append(mse)
            iteration_metrics['composite'].append(composite)
            
            print(f"Iteration {iter_idx}: Comp={composite:.6f} P={pearson:.4f} S={spearman:.4f} C={cosine:.4f} M={mse:.4f}")
        else:
            print(f"Warning: final_iteration_{iter_idx}.txt missing or empty")
            iteration_images.append(np.zeros(TOTAL, dtype=np.uint8))
            for metric in iteration_metrics:
                iteration_metrics[metric].append(0.0)
    
    # Calculate ideal case metrics
    ideal_metrics = None
    if ideal_image.size > 0:
        ideal_pearson, ideal_spearman, ideal_cosine, ideal_mse, ideal_composite = calculate_all_metrics(ideal_image, ref_image)
        ideal_metrics = {
            'pearson': ideal_pearson,
            'spearman': ideal_spearman,
            'cosine': ideal_cosine,
            'mse': ideal_mse,
            'composite': ideal_composite
        }
        print(f"Ideal case: Comp={ideal_composite:.6f} P={ideal_pearson:.4f} S={ideal_spearman:.4f} C={ideal_cosine:.4f} M={ideal_mse:.4f}")
    
    # Load summary data from Verilog
    summary_data = load_summary_data()
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    create_metric_comparison_plots(iteration_metrics, ideal_metrics)
    create_image_gallery(iteration_images, ref_image, ideal_image)
    
    # Generate detailed report
    print("Generating detailed report...")
    generate_detailed_report(iteration_metrics, ideal_metrics, summary_data)
    
    # Final summary
    print("\n" + "=" * 60)
    print("MULTI-METRIC ATTACK ANALYSIS COMPLETE")
    print("=" * 60)
    print("Generated files:")
    print("• multi_metric_analysis.png      - Comprehensive metric analysis")
    print("• multi_metric_gallery.png       - Visual results gallery")
    print("• multi_metric_detailed_report.txt - Complete analysis report")
    
    if len(iteration_metrics['composite']) > 0:
        best_idx = np.argmax(iteration_metrics['composite'])
        print(f"\nBEST RESULT: Iteration {best_idx} with {iteration_metrics['composite'][best_idx]:.6f} composite score")
        
        if ideal_metrics:
            effectiveness = (iteration_metrics['composite'][best_idx] / ideal_metrics['composite']) * 100
            print(f"EFFECTIVENESS: {effectiveness:.1f}% of theoretical maximum")
    
    print("\nMulti-metric analysis complete! 🎯")

if __name__ == "__main__":
    main()