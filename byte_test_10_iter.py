import warnings, glob, re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

W, H = 128, 128
TOTAL = W * H
NUM_ITERATIONS = 10

# ----------------------------------------------------------------------
#  Helper Functions
# ----------------------------------------------------------------------
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

def corr(a, b):
    """Calculate Pearson correlation, handling edge cases"""
    if len(a) == 0 or len(b) == 0 or len(a) != len(b): 
        return 0.0
    try:
        c, _ = pearsonr(a[:TOTAL], b[:TOTAL])
        return 0.0 if np.isnan(c) else c
    except:
        return 0.0

def calculate_metrics(img, ref):
    """Calculate correlation and pixel match percentage"""
    if len(img) == 0 or len(ref) == 0:
        return 0.0, 0.0
    correlation = corr(img, ref)
    pixel_match = np.mean(img == ref) * 100
    return correlation, pixel_match

# ----------------------------------------------------------------------
#  Load Reference Images
# ----------------------------------------------------------------------
print("Loading reference images...")
ref = load_img(Path("reference_correct_decode.txt"))
if ref.size == 0:
    print("ERROR: reference_correct_decode.txt missing or empty!")
    exit(1)

# Load ideal comparison case (correct TRNG_A, wrong TRNG_D)
ideal_case = load_img(Path("correct_trnga_wrong_trngd.txt"))
if ideal_case.size == 0:
    print("Warning: correct_trnga_wrong_trngd.txt missing!")
    ideal_case = None

print(f"Reference image loaded: {ref.shape}")

# ----------------------------------------------------------------------
#  Load Multi-Iteration Results
# ----------------------------------------------------------------------
print("Loading multi-iteration attack results...")
iteration_images = []
iteration_correlations = []
iteration_matches = []

for iter_idx in range(NUM_ITERATIONS):
    img_path = Path(f"final_iteration_{iter_idx}.txt")
    img = load_img(img_path)
    
    if img.size > 0:
        corr_val, match_val = calculate_metrics(img, ref)
        iteration_images.append(img)
        iteration_correlations.append(corr_val)
        iteration_matches.append(match_val)
        print(f"Iteration {iter_idx}: Correlation {corr_val:.6f}, Match {match_val:.2f}%")
    else:
        print(f"Warning: final_iteration_{iter_idx}.txt missing or empty")
        iteration_images.append(np.zeros(TOTAL, dtype=np.uint8))
        iteration_correlations.append(0.0)
        iteration_matches.append(0.0)

# Calculate ideal case metrics
if ideal_case is not None:
    ideal_corr, ideal_match = calculate_metrics(ideal_case, ref)
    print(f"Ideal case: Correlation {ideal_corr:.6f}, Match {ideal_match:.2f}%")
else:
    ideal_corr, ideal_match = 0.0, 0.0

# ----------------------------------------------------------------------
#  1. Iteration Progress Visualization
# ----------------------------------------------------------------------
plt.figure(figsize=(14, 6))

# Correlation progress
plt.subplot(1, 2, 1)
iterations = list(range(NUM_ITERATIONS))
corr_values = iteration_correlations

bars = plt.bar(iterations, corr_values, alpha=0.7, color='skyblue', edgecolor='navy')
if ideal_case is not None:
    plt.axhline(y=ideal_corr, color='red', linestyle='--', linewidth=2, 
                label=f'Ideal Case ({ideal_corr:.4f})')

plt.xlabel('Iteration Number')
plt.ylabel('Pearson Correlation')
plt.title('Attack Convergence Across Iterations')
plt.xticks(iterations)
plt.ylim(-1.0, 1.0)
plt.grid(True, alpha=0.3)
plt.legend()

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, corr_values)):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
             f'{val:.4f}', ha='center', fontsize=9, rotation=90)

# Pixel match progress
plt.subplot(1, 2, 2)
match_bars = plt.bar(iterations, iteration_matches, alpha=0.7, color='lightgreen', edgecolor='darkgreen')
if ideal_case is not None:
    plt.axhline(y=ideal_match, color='red', linestyle='--', linewidth=2, 
                label=f'Ideal Case ({ideal_match:.1f}%)')

plt.xlabel('Iteration Number')
plt.ylabel('Pixel Match (%)')
plt.title('Pixel-Level Accuracy Across Iterations')
plt.xticks(iterations)
plt.ylim(0, 100)
plt.grid(True, alpha=0.3)
plt.legend()

# Add value labels on bars
for i, (bar, val) in enumerate(zip(match_bars, iteration_matches)):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
             f'{val:.1f}%', ha='center', fontsize=9)

plt.tight_layout()
plt.savefig("iteration_progress.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
#  2. Multi-Iteration Image Gallery
# ----------------------------------------------------------------------
# Create comprehensive gallery: iterations + reference + ideal case
n_cols = 4
n_rows = 3
fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 12))
axes = axes.ravel()

# Show iteration results (first 10 slots)
for i in range(min(NUM_ITERATIONS, 10)):
    if i < len(iteration_images):
        axes[i].imshow(iteration_images[i].reshape(H, W), cmap="gray", vmin=0, vmax=255)
        axes[i].set_title(f"Iteration {i}\nCorr: {iteration_correlations[i]:.4f}\nMatch: {iteration_matches[i]:.1f}%", 
                         fontsize=10)
    else:
        axes[i].text(0.5, 0.5, "No Data", ha='center', va='center', transform=axes[i].transAxes)
    axes[i].axis("off")

# Show reference image
if len(axes) > 10:
    axes[10].imshow(ref.reshape(H, W), cmap="gray", vmin=0, vmax=255)
    axes[10].set_title("Reference\n(Original)", fontsize=10, color='black', weight='bold')
    axes[10].axis("off")

# Show ideal case
if len(axes) > 11 and ideal_case is not None:
    axes[11].imshow(ideal_case.reshape(H, W), cmap="gray", vmin=0, vmax=255)
    axes[11].set_title(f"Ideal Case\nCorr: {ideal_corr:.4f}\nMatch: {ideal_match:.1f}%", 
                     fontsize=10, color='red', weight='bold')
    axes[11].axis("off")
elif len(axes) > 11:
    axes[11].text(0.5, 0.5, "Ideal Case\nNot Available", ha='center', va='center', 
                 transform=axes[11].transAxes)
    axes[11].axis("off")

fig.suptitle("Multi-Iteration Attack Results Gallery", weight="bold", fontsize=14)
fig.tight_layout()
plt.savefig("multi_iteration_gallery.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
#  3. Detailed Comparison: Best vs Worst vs Ideal
# ----------------------------------------------------------------------
if len(iteration_correlations) > 0:
    best_idx = np.argmax(iteration_correlations)
    worst_idx = np.argmin(iteration_correlations)
    
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    
    # Original reference
    axes[0].imshow(ref.reshape(H, W), cmap="gray", vmin=0, vmax=255)
    axes[0].set_title("Original\nReference", fontsize=12, weight='bold')
    axes[0].axis("off")
    
    # Worst iteration
    axes[1].imshow(iteration_images[worst_idx].reshape(H, W), cmap="gray", vmin=0, vmax=255)
    axes[1].set_title(f"Worst Result\n(Iter {worst_idx})\nCorr: {iteration_correlations[worst_idx]:.4f}", 
                     fontsize=12, color='red')
    axes[1].axis("off")
    
    # Best iteration
    axes[2].imshow(iteration_images[best_idx].reshape(H, W), cmap="gray", vmin=0, vmax=255)
    axes[2].set_title(f"Best Result\n(Iter {best_idx})\nCorr: {iteration_correlations[best_idx]:.4f}", 
                     fontsize=12, color='green')
    axes[2].axis("off")
    
    # Ideal case
    if ideal_case is not None:
        axes[3].imshow(ideal_case.reshape(H, W), cmap="gray", vmin=0, vmax=255)
        axes[3].set_title(f"Ideal Case\n(Correct TRNG_A)\nCorr: {ideal_corr:.4f}", 
                         fontsize=12, color='blue')
        axes[3].axis("off")
    else:
        axes[3].text(0.5, 0.5, "Ideal Case\nNot Available", ha='center', va='center', 
                    transform=axes[3].transAxes)
        axes[3].axis("off")
    
    # Difference map (best vs ideal)
    if ideal_case is not None:
        diff = np.abs(iteration_images[best_idx].astype(int) - ideal_case.astype(int))
        axes[4].imshow(diff.reshape(H, W), cmap="hot", vmin=0, vmax=255)
        axes[4].set_title(f"Difference Map\n(Best vs Ideal)\nMean: {np.mean(diff):.1f}", fontsize=12)
        axes[4].axis("off")
    else:
        axes[4].text(0.5, 0.5, "Difference\nNot Available", ha='center', va='center', 
                    transform=axes[4].transAxes)
        axes[4].axis("off")
    
    fig.suptitle("Attack Quality Analysis: Best vs Worst vs Ideal", weight="bold", fontsize=14)
    fig.tight_layout()
    plt.savefig("attack_quality_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

# ----------------------------------------------------------------------
#  4. Convergence Analysis
# ----------------------------------------------------------------------
plt.figure(figsize=(12, 8))

# Plot 1: Correlation trend
plt.subplot(2, 2, 1)
plt.plot(range(NUM_ITERATIONS), iteration_correlations, 
         'o-', linewidth=2, markersize=8, color='blue')
if ideal_case is not None:
    plt.axhline(y=ideal_corr, color='red', linestyle='--', 
                label=f'Theoretical Max ({ideal_corr:.4f})')
plt.xlabel('Iteration')
plt.ylabel('Correlation')
plt.title('Correlation Convergence')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xticks(range(NUM_ITERATIONS))
plt.ylim(-1.0, 1.0)

# Plot 2: Pixel match trend
plt.subplot(2, 2, 2)
plt.plot(range(NUM_ITERATIONS), iteration_matches, 
         'o-', linewidth=2, markersize=8, color='green')
if ideal_case is not None:
    plt.axhline(y=ideal_match, color='red', linestyle='--', 
                label=f'Theoretical Max ({ideal_match:.1f}%)')
plt.xlabel('Iteration')
plt.ylabel('Pixel Match (%)')
plt.title('Pixel Match Convergence')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xticks(range(NUM_ITERATIONS))

# Plot 3: Improvement analysis
plt.subplot(2, 2, 3)
if len(iteration_correlations) > 1:
    improvements = [iteration_correlations[i+1] - iteration_correlations[i] 
                   for i in range(len(iteration_correlations)-1)]
    plt.bar(range(1, NUM_ITERATIONS), improvements, alpha=0.7, 
            color=['green' if imp > 0 else 'red' for imp in improvements])
    plt.xlabel('Iteration Transition')
    plt.ylabel('Correlation Improvement')
    plt.title('Iteration-to-Iteration Improvement')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(1, NUM_ITERATIONS), [f'{i-1}→{i}' for i in range(1, NUM_ITERATIONS)])

# Plot 4: Attack effectiveness
plt.subplot(2, 2, 4)
if ideal_case is not None and ideal_corr > 0:
    effectiveness = [(c/ideal_corr)*100 for c in iteration_correlations]
    plt.plot(range(NUM_ITERATIONS), effectiveness, 'o-', linewidth=2, markersize=8, color='purple')
    plt.xlabel('Iteration')
    plt.ylabel('Attack Effectiveness (%)')
    plt.title('Effectiveness vs Theoretical Maximum')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(NUM_ITERATIONS))
    plt.ylim(0, 100)

plt.tight_layout()
plt.savefig("convergence_analysis.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
#  5. Statistical Summary Report
# ----------------------------------------------------------------------
with open("multi_iteration_report.txt", "w") as fh:
    fh.write("Multi-Iteration Cryptographic Attack Analysis\n")
    fh.write("=" * 50 + "\n\n")
    
    # Attack overview
    fh.write("ATTACK OVERVIEW:\n")
    fh.write("-" * 16 + "\n")
    fh.write(f"Number of iterations: {NUM_ITERATIONS}\n")
    fh.write(f"Image size: {W}×{H} = {TOTAL} pixels\n\n")
    
    # Iteration results
    fh.write("ITERATION RESULTS:\n")
    fh.write("-" * 18 + "\n")
    fh.write("Iter |  Correlation  |  Pixel Match  |  Status\n")
    fh.write("-----|---------------|---------------|----------\n")
    
    for i in range(NUM_ITERATIONS):
        if i < len(iteration_correlations):
            status = "BEST" if i == np.argmax(iteration_correlations) else ""
            status = "WORST" if i == np.argmin(iteration_correlations) else status
            fh.write(f" {i:2d}  |    {iteration_correlations[i]:8.6f}    |    {iteration_matches[i]:6.2f}%    | {status:>8s}\n")
        else:
            fh.write(f" {i:2d}  |     MISSING     |     MISSING     |   ERROR\n")
    
    fh.write("\n")
    
    # Statistical analysis
    if len(iteration_correlations) > 0:
        fh.write("STATISTICAL ANALYSIS:\n")
        fh.write("-" * 21 + "\n")
        fh.write(f"Best correlation:     {max(iteration_correlations):8.6f} (Iteration {np.argmax(iteration_correlations)})\n")
        fh.write(f"Worst correlation:    {min(iteration_correlations):8.6f} (Iteration {np.argmin(iteration_correlations)})\n")
        fh.write(f"Mean correlation:     {np.mean(iteration_correlations):8.6f}\n")
        fh.write(f"Std deviation:        {np.std(iteration_correlations):8.6f}\n")
        
        if len(iteration_correlations) > 1:
            improvement = iteration_correlations[-1] - iteration_correlations[0]
            fh.write(f"Total improvement:    {improvement:+8.6f}\n")
        
        fh.write("\n")
    
    # Comparison with ideal case
    if ideal_case is not None:
        fh.write("COMPARISON WITH IDEAL CASE:\n")
        fh.write("-" * 28 + "\n")
        fh.write(f"Ideal correlation:    {ideal_corr:8.6f}\n")
        fh.write(f"Ideal pixel match:    {ideal_match:6.2f}%\n")
        
        if len(iteration_correlations) > 0:
            best_corr = max(iteration_correlations)
            effectiveness = (best_corr / ideal_corr) * 100 if ideal_corr > 0 else 0
            fh.write(f"Attack effectiveness: {effectiveness:6.2f}% of theoretical maximum\n")
        
        fh.write("\n")
    
    # Attack convergence
    if len(iteration_correlations) > 1:
        fh.write("CONVERGENCE ANALYSIS:\n")
        fh.write("-" * 21 + "\n")
        
        # Check if attack is converging
        last_half = iteration_correlations[NUM_ITERATIONS//2:]
        if len(last_half) > 1:
            trend = np.polyfit(range(len(last_half)), last_half, 1)[0]
            if trend > 0.001:
                convergence = "IMPROVING"
            elif trend < -0.001:
                convergence = "DEGRADING"
            else:
                convergence = "STABLE"
            fh.write(f"Convergence trend:    {convergence}\n")
        
        # Find best improvement step
        improvements = [iteration_correlations[i+1] - iteration_correlations[i] 
                       for i in range(len(iteration_correlations)-1)]
        best_step = np.argmax(improvements)
        fh.write(f"Best improvement:     Iteration {best_step} → {best_step+1} ({improvements[best_step]:+.6f})\n")
        
        fh.write("\n")
    
    # Recommendations
    fh.write("RECOMMENDATIONS:\n")
    fh.write("-" * 16 + "\n")
    
    if len(iteration_correlations) > 0:
        best_corr = max(iteration_correlations)
        if best_corr < 0.5:
            fh.write("• Attack showing poor results - consider different approach\n")
        elif best_corr < 0.8:
            fh.write("• Moderate success - may need more iterations or refinement\n")
        else:
            fh.write("• High success rate - attack is effective\n")
        
        if len(iteration_correlations) > 1:
            if iteration_correlations[-1] > iteration_correlations[0]:
                fh.write("• Iterations are helping - continue refinement\n")
            else:
                fh.write("• No improvement from iterations - check algorithm\n")
    
    if ideal_case is not None and len(iteration_correlations) > 0:
        effectiveness = (max(iteration_correlations) / ideal_corr) * 100 if ideal_corr > 0 else 0
        if effectiveness < 50:
            fh.write("• Large gap vs theoretical maximum - algorithm needs improvement\n")
        elif effectiveness < 80:
            fh.write("• Reasonable performance but room for optimization\n")
        else:
            fh.write("• Excellent performance - close to theoretical maximum\n")

# ----------------------------------------------------------------------
#  Summary Output
# ----------------------------------------------------------------------
print("\n" + "="*60)
print("MULTI-ITERATION ATTACK ANALYSIS COMPLETE")
print("="*60)
print("Generated files:")
print("• iteration_progress.png        - Attack convergence charts")
print("• multi_iteration_gallery.png   - Complete image gallery") 
print("• attack_quality_comparison.png - Best vs worst vs ideal")
print("• convergence_analysis.png      - Detailed convergence analysis")
print("• multi_iteration_report.txt    - Comprehensive text report")

if len(iteration_correlations) > 0:
    best_idx = np.argmax(iteration_correlations)
    print(f"\nBEST RESULT: Iteration {best_idx} with {iteration_correlations[best_idx]:.6f} correlation")
    
    if ideal_case is not None:
        effectiveness = (iteration_correlations[best_idx] / ideal_corr) * 100 if ideal_corr > 0 else 0
        print(f"EFFECTIVENESS: {effectiveness:.1f}% of theoretical maximum")

print("\nAnalysis complete! 🎯")