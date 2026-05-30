import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('/media/jannik/Daten/Dokumente/UNI/Doktor/ContinuousDijkstra/results/boxplot_data.csv', sep="\t")
print(list(df.columns))
print(df.head())

pa = [df.loc[df['type'] == 'PA']['time_ms'].values]
fastes = [df.loc[df['type'] == 'Fastes']['time_ms'].values]
noDeadend = [df.loc[df['type'] == 'FastesNoDeadend']['time_ms'].values]
print(np.mean(pa), np.mean(fastes), np.mean(noDeadend))
print(np.median(pa), np.median(fastes), np.median(noDeadend))
print(np.std(pa), np.std(fastes), np.std(noDeadend))
print(np.min(pa), np.min(fastes), np.min(noDeadend))
print(np.max(pa), np.max(fastes), np.max(noDeadend))

# Create boxplot
plt.figure(figsize=(4, 6))
plt.boxplot([df.loc[df['type'] == 'Fastes']['time_ms'].values],
             labels=['Fastes'])
plt.boxplot([df.loc[df['type'] == 'FastesNoDeadend']['time_ms'].values], 
             labels=['FastesNoDeadend'])
plt.boxplot([df.loc[df['type'] == 'PA']['time_ms'].values], 
             labels=[ 'PA'])
plt.ylabel('Values')
plt.xlabel('Groups')
plt.title('Boxplot Example')
plt.grid(axis='y', alpha=0.3)
plt.show()