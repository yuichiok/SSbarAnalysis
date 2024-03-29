from ROOT import TH2F
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go


particles = ["Pi", "K", "p", "mu", "e"]
# Read TH2F histogram from ROOT file
file = uproot.open("dedx2violin.root")

# Extract x, y, and z values as numpy arrays
data = {}
dfc  = {}
for particle in particles:
  h = file[f"h2_dd_reco_Pi_{particle}_nocut_dEdx_p"]
  data[particle] = {
    "x": h.axis(0).edges()[:-1],
    "y": h.axis(1).edges()[:-1],
    "z": h.values().flatten(),
  }
  x_values = data[particle]["x"]
  y_values = data[particle]["y"]
  z_values = data[particle]["z"]

  dfc[f"{particle}_x"] = pd.Series(np.repeat(x_values, len(y_values)))
  dfc[f"{particle}_y"] = pd.Series(np.tile(y_values, len(x_values)))
  dfc[f"{particle}_z"] = pd.Series(z_values)

  # dfs[particle] = pd.DataFrame({
  #   "x": pd.Series(np.repeat(x_values, len(y_values))),
  #   "y": pd.Series(np.tile(y_values, len(x_values))),
  #   "z": pd.Series(z_values)
  # })

df = pd.DataFrame(dfc)

fig = go.Figure()

for particle in particles:
  fig.add_trace(go.Violin(x=df[f"{particle}_x"],
                          y=df[f"{particle}_y"],
                          legendgroup='Yes', scalegroup='Yes', name='Yes',
                          side='positive',
                          line_color='blue')
               )


# fig.add_trace(go.Violin(x=df['x'],
#                         y=df['y'],
#                         legendgroup='Yes', scalegroup='Yes', name='Yes',
#                         side='negative',
#                         line_color='blue')
#              )

fig.update_traces(meanline_visible=True)
fig.update_layout(violingap=0, violinmode='overlay')
fig.show()

# Display the DataFrame
print(df)
