import re
import os
import sys
import matplotlib.pyplot as plt

def extract_sample_name_and_assigned_rate(log_file):
    """
    Extracts sample name and assigned rate from the given log file.

    Args:
        log_file: Path to the log file.

    Returns:
        tuple: A tuple containing the sample name and assigned rate.
    """
    with open(log_file, 'r') as f:
        log_content = f.read()

    # Extract sample name
    output_file_match = re.search(r"Output file : (\S+)\.fC", log_content)
    if output_file_match:
        sample_name = output_file_match.group(1)
    else:
        sample_name = "Unknown"

    # Extract assigned rate
    assigned_rate_match = re.search(r"Successfully assigned alignments : (\d+) \((\d+\.\d+)%\)", log_content)
    if assigned_rate_match:
        assigned_rate = float(assigned_rate_match.group(2))
    else:
        assigned_rate = None

    return sample_name, assigned_rate

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <log_path>")
        sys.exit(1)

    log_path = sys.argv[1]

    sample_rates = []
    sample_names = []

    for filename in os.listdir(log_path):
        if filename.startswith("fe"):
            log_file = os.path.join(log_path, filename)
            sample_name, assigned_rate = extract_sample_name_and_assigned_rate(log_file)

            if assigned_rate is not None:
                sample_rates.append(assigned_rate)
                sample_names.append(sample_name)

    # Sort samples by assigned rate
    sorted_indices = sorted(range(len(sample_rates)), key=lambda k: sample_rates[k])
    sorted_sample_rates = [sample_rates[i] for i in sorted_indices]
    sorted_sample_names = [sample_names[i] for i in sorted_indices]

    # Write sorted results to file
    with open("summary.txt", "w") as f:
        f.write("Sample\tAssigned Rate\n")
        for sample, rate in zip(sorted_sample_names, sorted_sample_rates):
            f.write(f"{sample}\t{rate:.2f}%\n")

    # Create box plot
    plt.figure(figsize=(5, 6))
    plt.boxplot(sorted_sample_rates, showfliers=True)
    plt.xlabel("Samples")
    plt.ylabel("Assigned Rate (%)")
    plt.title("Assigned Rate Distribution")

    # Add labels for 5 lowest samples
    hapos = ['left','right']
    for i in range(5):
        plt.text(1, sorted_sample_rates[i], sorted_sample_names[i], ha=hapos[i%2], va='bottom')

    plt.savefig("assigned_rate_boxplot.png")
    plt.show()

if __name__ == "__main__":
    main()
