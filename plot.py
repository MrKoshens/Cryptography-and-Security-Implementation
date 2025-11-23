import matplotlib.pyplot as plt
import seaborn as sns
import re

sns.set(style="whitegrid")

SORT_NAMES = ["Bubble", "Merge", "Quick", "Heap"]
METRICS   = ["Swaps", "Comps", "Cycles"]
STAT_TYPES = ["Min", "Max", "Median", "Avg"]
SIZES = list(range(100, 1001, 100))

# stats[metric][stat_type][sort_name] = list of values
stats = {
    metric: {
        stat: {name: [] for name in SORT_NAMES}
        for stat in STAT_TYPES
    }
    for metric in METRICS
}

# Read entire file into memory
with open("results.txt", "r") as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith("Array Size:"):
        i += 1
        for sort in SORT_NAMES:
            header = lines[i].strip()
            assert header.startswith(f"{sort} Sort:"), f"Expected '{sort} Sort:' but got '{header}'"
            i += 1
            for metric in METRICS:
                stat_line = lines[i].strip()
                m = re.search(
                    r"Min=([\d\.]+)\s+Max=([\d\.]+)\s+Median=([\d\.]+)\s+Avg=([\d\.]+)",
                    stat_line
                )
                if not m:
                    raise ValueError(f"Could not parse stats from: '{stat_line}'")
                min_v, max_v, med_v, avg_v = map(float, m.groups())
                stats[metric]["Min"][sort].append(min_v)
                stats[metric]["Max"][sort].append(max_v)
                stats[metric]["Median"][sort].append(med_v)
                stats[metric]["Avg"][sort].append(avg_v)
                i += 1
            # skip blank line
            i += 1
    else:
        i += 1

def plot_stat(metric, stat_type):
    plt.figure(figsize=(8, 5))
    for sort in SORT_NAMES:
        plt.plot(
            SIZES,
            stats[metric][stat_type][sort],
            marker="o",
            label=sort
        )
    plt.title(f"{metric} — {stat_type} vs Array Size")
    plt.xlabel("Array Size")
    plt.ylabel(f"{metric} ({stat_type})")
    plt.legend()
    plt.tight_layout()
    filename = f"{metric.lower()}_{stat_type.lower()}.png"
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {filename}")

for metric in METRICS:
    for stat_type in STAT_TYPES:
        plot_stat(metric, stat_type)