import plotly.express as px
import pandas as pd
dpPTM = pd.read_csv("../data/dpPTMextended.tsv", sep="\t")
x = dpPTM["conservation"].to_numpy()
y = dpPTM["sequonEfficiency"].to_numpy()
print(x)
print(y)
fig = px.scatter(x=x, y=y, trendline="ols")
fig.show()