# To test how "df.loc[len(df)] = [proc, d.sumEntries(), 1, 1]" work

import pandas as pd

# Initialize a DataFrame with the specified columns
df = pd.DataFrame(columns=['proc', 'sumEntries', 'nRV', 'nWV'])

# Simulated process list and dataset (in reality, these would come from RooFit)
processes = ['proc1', 'proc2', 'proc3']

# Example of a class simulating what d.sumEntries() would do
class FakeDataset:
    def __init__(self, entries):
        self.entries = entries
    
    def sumEntries(self):
        # Simulate the sumEntries method
        return self.entries

# Simulated dataset results for each process
datasets = {
    'proc1': FakeDataset(100),
    'proc2': FakeDataset(250),
    'proc3': FakeDataset(400)
}

# Loop over processes and add entries to the DataFrame
for proc in processes:
    d = datasets[proc]  # Here, d is the simulated dataset for the current process
    df.loc[len(df)] = [proc, d.sumEntries(), 1, 1]

# Output the DataFrame
print(df)
