import sys, re
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

def parse_fastp_log(log_string):
    """
    Parses a fastp log string and extracts relevant information.

    Args:
        log_string: The content of the fastp log file as a string.

    Returns:
        A dictionary containing the extracted information.
    """
    data = {}
    # Extract Sample_name
    match = re.search(r"HTML report: (.*).fastp.html", log_string)
    if match:
        data['Sample_name'] = match.group(1)
    # Extract Raw_reads
    match = re.search(r"Read1 before filtering:\n.*?total reads: (\d+)", log_string, re.DOTALL)
    if match:
        data['Raw_reads'] = int(match.group(1))

    # Extract Clean_reads and Q20%, Q30%
    match = re.search(r"Read1 after filtering:\n.*?total reads: (\d+)\n.*?Q20 bases: \d+\((\d+\.\d+)\%\)\n.*?Q30 bases: \d+\((\d+\.\d+)\%\)", log_string, re.DOTALL)
    if match:
        data['Clean_reads'] = int(match.group(1))
        data['Q20%'] = float(match.group(2))
        data['Q30%'] = float(match.group(3))

    # Extract Duplication_rate
    match = re.search(r"Duplication rate: (\d+\.\d+)\%", log_string)
    if match:
        data['Duplication_rate'] = float(match.group(1))
    # Extract Insert peak
    match = re.search(r"Insert size peak \(evaluated by paired-end reads\): (\d+)", log_string)
    if match:
        data['Insert_peak'] = float(match.group(1))
    return data

def summarize_fastp_logs(input_path):
    """
    Summarizes multiple fastp log files in a given directory.

    Args:
        input_path: Path to the directory containing the log files.

    Returns:
        A pandas DataFrame containing the summarized data.
    """
    log_files = list(Path(input_path).glob('fastp.*'))
    data = []
    for log_file in log_files:
        print(log_file)
        with open(log_file, 'r') as f:
            log_string = f.read()
        data.append(parse_fastp_log(log_string))
    df = pd.DataFrame(data)
    df_sorted_asc = df.sort_values(by='Sample_name')
    return df_sorted_asc

def plot_duplication_rates(df):
    """
    Plots a vertical boxplot of duplication rates for all samples and annotates the top three highest points.

    Args:
        df: A pandas DataFrame containing the summarized data.
    """
    plt.figure(figsize=(4, 5))
    plt.boxplot(df['Duplication_rate'], vert=True)
    plt.title('Duplication Rates Across Samples')
    plt.ylabel('Duplication Rate (%)')

    # Annotate the top three highest points
    y = df['Duplication_rate']
    x = [1] * len(y)  # All points are in the first boxplot
    for i, v in enumerate(sorted(zip(y, df['Sample_name']), reverse=True)[:3]):
        plt.text(x[i], v[0], f"{v[1]}: {v[0]:.2f}%", ha='right', va='bottom')

    plt.savefig('duplication_rates_boxplot.png')
    plt.show()


input_dir = sys.argv[1]
summary_df = summarize_fastp_logs(input_dir)

# Write summary to a text file
with open('sequencing_summary.txt', 'w') as f:
    f.write(summary_df.to_string(index=False))

print("Summary written to fastp_summary.txt")

# Plot duplication rates and save the plot
plot_duplication_rates(summary_df)
