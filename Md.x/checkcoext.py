import numpy as np
from scipy import signal
from scipy.stats import kurtosis, skew
import matplotlib.pyplot as plt

def analyze_coexistence(density_file='density.dat', plot=False):
    """
    Comprehensive solid-liquid coexistence analysis using adaptive thresholds
    
    Parameters:
    -----------
    density_file : str
        Path to density profile file
    plot : bool
        Whether to generate visualization plots
    
    Returns:
    --------
    dict : Analysis results and metrics
    """
    
    # Load and validate data
    try:
        data = np.loadtxt(density_file)
        positions = data[:, 0]
        densities = data[:, 1]
    except:
        raise FileNotFoundError(f"Cannot read {density_file}")
    
    # Basic statistics
    mean_density = np.mean(densities)
    std_density = np.std(densities)
    cv = std_density / mean_density  # Coefficient of variation
    
    # Advanced statistical measures
    density_skewness = skew(densities)
    density_kurtosis = kurtosis(densities)
    bimodality_coeff = (density_skewness**2 + 1) / (density_kurtosis + 3)
    
    # Adaptive thresholds using multiple methods
    # Method 1: Statistical (Otsu-like)
    k_factor = 1.5
    high_threshold = mean_density + k_factor * std_density
    low_threshold = mean_density - k_factor * std_density
    
    # Method 2: Percentile-based thresholds
    p25, p75 = np.percentile(densities, [25, 75])
    iqr = p75 - p25
    outlier_high = p75 + 1.5 * iqr
    outlier_low = p25 - 1.5 * iqr
    
    # Use the more conservative threshold
    final_high = min(high_threshold, outlier_high)
    final_low = max(low_threshold, outlier_low)
    
    # Phase identification
    solid_mask = densities > final_high
    liquid_mask = densities < final_low
    interface_mask = ~(solid_mask | liquid_mask)
    
    # Calculate phase fractions
    solid_fraction = np.sum(solid_mask) / len(densities)
    liquid_fraction = np.sum(liquid_mask) / len(densities)
    interface_fraction = np.sum(interface_mask) / len(densities)
    
    # Peak detection for bimodality
    hist, bin_edges = np.histogram(densities, bins=max(20, len(densities)//20))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Find peaks in histogram
    peaks, properties = signal.find_peaks(hist, 
                                        height=len(densities)*0.02,  # At least 2% of data
                                        distance=len(hist)//10)       # Minimum separation
    
    n_peaks = len(peaks)
    bimodal = n_peaks >= 2
    
    # Interface analysis
    gradients = np.abs(np.gradient(densities))
    mean_gradient = np.mean(gradients)
    max_gradient = np.max(gradients)
    
    # Sharp interface detection
    sharp_interfaces = np.where(gradients > mean_density * 0.1)[0]
    n_interfaces = len(sharp_interfaces)
    
    # Coexistence scoring system (0-6 scale)
    scores = {
        'bimodal': int(bimodal),
        'high_cv': int(cv > 0.5),
        'phase_balance': int(solid_fraction > 0.1 and liquid_fraction > 0.1),
        'sharp_gradients': int(mean_gradient > std_density),
        'bimodality_coeff': int(bimodality_coeff > 0.55),
        'multiple_interfaces': int(n_interfaces > 2)
    }
    
    coexistence_score = sum(scores.values())
    
    # System state classification
    if coexistence_score >= 4:
        system_state = 'STRONG SOLID-LIQUID COEXISTENCE'
    elif coexistence_score >= 3:
        system_state = 'SOLID-LIQUID COEXISTENCE'
    elif coexistence_score >= 2:
        system_state = 'POSSIBLE COEXISTENCE'
    elif solid_fraction > 0.8:
        system_state = 'PREDOMINANTLY SOLID'
    elif liquid_fraction > 0.8:
        system_state = 'PREDOMINANTLY LIQUID'
    else:
        system_state = 'SINGLE PHASE'
    
    # Phase labels
    phase_labels = np.full(len(densities), 'INTERFACE', dtype=object)
    phase_labels[solid_mask] = 'SOLID'
    phase_labels[liquid_mask] = 'LIQUID'
    
    # Continuous regions analysis
    solid_regions = find_continuous_regions(solid_mask)
    liquid_regions = find_continuous_regions(liquid_mask)
    
    results = {
        'system_state': system_state,
        'coexistence_score': coexistence_score,
        'score_breakdown': scores,
        'statistics': {
            'mean_density': mean_density,
            'std_density': std_density,
            'cv': cv,
            'skewness': density_skewness,
            'kurtosis': density_kurtosis,
            'bimodality_coefficient': bimodality_coeff
        },
        'thresholds': {
            'high_threshold': final_high,
            'low_threshold': final_low
        },
        'fractions': {
            'solid': solid_fraction,
            'liquid': liquid_fraction,
            'interface': interface_fraction
        },
        'interface_analysis': {
            'n_peaks': n_peaks,
            'mean_gradient': mean_gradient,
            'max_gradient': max_gradient,
            'n_interfaces': n_interfaces
        },
        'regions': {
            'solid_regions': solid_regions,
            'liquid_regions': liquid_regions
        },
        'data': {
            'positions': positions,
            'densities': densities,
            'phase_labels': phase_labels
        }
    }
    
    # Print summary
    print_summary(results)
    
    # Optional plotting
    if plot:
        create_plots(results)
    
    return results

def find_continuous_regions(mask):
    """Find continuous regions where mask is True"""
    regions = []
    start = None
    
    for i, val in enumerate(mask):
        if val and start is None:
            start = i
        elif not val and start is not None:
            regions.append((start, i-1))
            start = None
    
    if start is not None:
        regions.append((start, len(mask)-1))
    
    return regions

def print_summary(results):
    """Print formatted analysis summary"""
    print("=" * 50)
    print("COEXISTENCE ANALYSIS RESULTS")
    print("=" * 50)
    print(f"System State: {results['system_state']}")
    print(f"Coexistence Score: {results['coexistence_score']}/6")
    print()
    
    print("STATISTICAL PARAMETERS:")
    stats = results['statistics']
    print(f"  Mean Density: {stats['mean_density']:.3f}")
    print(f"  Std Deviation: {stats['std_density']:.3f}")
    print(f"  Coeff. of Variation: {stats['cv']:.3f}")
    print(f"  Bimodality Coefficient: {stats['bimodality_coefficient']:.3f}")
    print()
    
    print("PHASE FRACTIONS:")
    fractions = results['fractions']
    print(f"  Solid: {fractions['solid']:.3f} ({fractions['solid']*100:.1f}%)")
    print(f"  Liquid: {fractions['liquid']:.3f} ({fractions['liquid']*100:.1f}%)")
    print(f"  Interface: {fractions['interface']:.3f} ({fractions['interface']*100:.1f}%)")
    print()
    
    print("INTERFACE ANALYSIS:")
    interface = results['interface_analysis']
    print(f"  Density Peaks: {interface['n_peaks']}")
    print(f"  Sharp Interfaces: {interface['n_interfaces']}")
    print(f"  Mean Gradient: {interface['mean_gradient']:.3f}")

def create_plots(results):
    """Create visualization plots"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    positions = results['data']['positions']
    densities = results['data']['densities']
    phase_labels = results['data']['phase_labels']
    
    # Plot 1: Density profile with phase regions
    colors = {'SOLID': 'red', 'LIQUID': 'blue', 'INTERFACE': 'gray'}
    for phase in ['SOLID', 'LIQUID', 'INTERFACE']:
        mask = phase_labels == phase
        ax1.scatter(positions[mask], densities[mask], 
                   c=colors[phase], label=phase, alpha=0.7, s=10)
    
    ax1.axhline(results['thresholds']['high_threshold'], 
                color='red', linestyle='--', alpha=0.5)
    ax1.axhline(results['thresholds']['low_threshold'], 
                color='blue', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Density')
    ax1.set_title('Density Profile with Phase Classification')
    ax1.legend()
    
    # Plot 2: Density histogram
    ax2.hist(densities, bins=30, alpha=0.7, edgecolor='black')
    ax2.axvline(results['thresholds']['high_threshold'], 
                color='red', linestyle='--', label='High threshold')
    ax2.axvline(results['thresholds']['low_threshold'], 
                color='blue', linestyle='--', label='Low threshold')
    ax2.set_xlabel('Density')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Density Distribution')
    ax2.legend()
    
    # Plot 3: Gradient analysis
    gradients = np.abs(np.gradient(densities))
    ax3.plot(positions, gradients, 'g-', alpha=0.7)
    ax3.axhline(results['interface_analysis']['mean_gradient'], 
                color='orange', linestyle='--', label='Mean gradient')
    ax3.set_xlabel('Position')
    ax3.set_ylabel('|Gradient|')
    ax3.set_title('Interface Sharpness Analysis')
    ax3.legend()
    
    # Plot 4: Score breakdown
    scores = results['score_breakdown']
    ax4.bar(scores.keys(), scores.values())
    ax4.set_ylabel('Score')
    ax4.set_title('Coexistence Criteria Breakdown')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()

# Usage examples:
if __name__ == "__main__":
    # Basic analysis
    results = analyze_coexistence('density.dat')
    
    # With plotting
    results = analyze_coexistence('density.dat', plot=True)
    
    # Access specific results
    print(f"\nSystem classification: {results['system_state']}")
    print(f"Confidence score: {results['coexistence_score']}/6")
