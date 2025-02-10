import os
import sys
import re
import matplotlib.pyplot as plt

def extract_info_from_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    
    sample_name = os.path.basename(os.path.dirname(file_path))
    
    total_reads_match = re.search(r'(\d+) reads; of these:', content)
    exact_match_rate_match = re.search(r'\d+ \((\d+\.\d+)\%\) aligned concordantly exactly 1 time', content)
    multiple_match_rate_match = re.search(r'\d+ \((\d+\.\d+)\%\) aligned concordantly >1 times', content)
    
    total_reads = int(total_reads_match.group(1)) if total_reads_match else None
    exact_match_rate = float(exact_match_rate_match.group(1).replace('%', '')) / 100 if exact_match_rate_match else None
    multiple_match_rate = float(multiple_match_rate_match.group(1).replace('%', '')) / 100 if multiple_match_rate_match else None
    
    return {
        'sample_name': sample_name,
        'total_reads': total_reads,
        'exact_match_rate': exact_match_rate,
        'multiple_match_rate': multiple_match_rate
    }

def main(input_path):
    summary_data = []
    
    # Sort directories by name
    dirs = sorted([d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d))])
    
    for dir_name in dirs:
        file_path = os.path.join(input_path, dir_name, f"{dir_name}.summary")
        if os.path.exists(file_path):
            info = extract_info_from_file(file_path)
            summary_data.append(info)
    
    # Write summary to text file
    with open('rna_seq_summary.txt', 'w') as f:
        f.write("Sample Name\tTotal Reads\tExact Match Rate\tMultiple Match Rate\n")
        for data in summary_data:
            f.write(f"{data['sample_name']}\t{data['total_reads']}\t{data['exact_match_rate']:.2%}\t{data['multiple_match_rate']:.2%}\n")
    
    # Plot exact match rates
    exact_match_rates = [d['exact_match_rate'] for d in summary_data]
    sample_names = [d['sample_name'] for d in summary_data]
    
    plt.figure(figsize=(4, 5))
    boxplot = plt.boxplot(exact_match_rates, vert=True, patch_artist=True, boxprops=dict(facecolor='blue'))
    
    # Get the positions of the whiskers and fliers
    whisker_positions = [whisker.get_ydata() for whisker in boxplot['whiskers']]
    flier_positions = [flier.get_ydata() for flier in boxplot['fliers']]
    
    # Combine whisker and flier positions to find the lowest points
    all_positions = whisker_positions + flier_positions
    
    # Find the five lowest exact match rates
    lowest_indices = sorted(range(len(exact_match_rates)), key=lambda k: exact_match_rates[k])[:5]
    
    # Add sample names near the five lowest points
    for idx in lowest_indices:
        plt.text(1.05, exact_match_rates[idx], sample_names[idx], fontsize=9, color='red', ha='left')
    plt.ylabel('Exact Match Rate')
    plt.title('RNA-Seq Exact Match Rates')
    plt.tight_layout()
    plt.savefig('exact_match_rates.png')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_data>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    main(input_path)
