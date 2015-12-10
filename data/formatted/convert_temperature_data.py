"""
Description
------------

Script to converts TempExp.csv file (the raw data file) into a different
format so it can be used in the IPM analysis.

Author
------
Mark Wilber

"""
import numpy as np
import pandas as pd

temp_data = pd.read_csv("../archival/TempExp.csv")

# Drop the cultures...though we could keep these if we want to compare frog
# to culture...also drop tadpoles and subadults at least for now

temp_trun = temp_data[temp_data.individual.str.contains('F.*')]

# 1. Split by individual

individuals = temp_trun.individual.unique()

ind_dfs = []

for indiv in individuals:

    tdat = temp_trun[temp_trun.individual == indiv]

    time_diff = np.array(tdat.day[1:]) - np.array(tdat.day[0:-1])

    matched_ZE = zip(tdat.ZE[0:-1], tdat.ZE[1:], time_diff,
                            np.repeat(indiv, len(time_diff)),
                            tdat.Temperature[:-1],
                            np.repeat(1, len(time_diff)), tdat.InitialStage[:-1])

    tdf = pd.DataFrame(matched_ZE, columns=['size', 'sizeNext', 'time',
                                'individual', 'temp', "surv", "stage"])

    ind_dfs.append(tdf)

full_data = pd.concat(ind_dfs)

# Drop the points where they are dead in both time points
full_data = full_data[~(np.array(full_data['size'] == 'dead'))]

# Change Dead to Nan
ind = full_data.sizeNext == "dead"
full_data.surv[ind] = 0
full_data.sizeNext[ind] = np.nan

# Change types to floats
full_data[['size', 'sizeNext']] = \
                full_data[['size', 'sizeNext']].astype(float)

# There is a weird aspect to the high temperature data where frogs would ramp
# up there loads and then reduce the load. For temp == 26, set the highest load
# to the load they die at

full_data.to_csv("converted_temp_data.csv", index=False)



