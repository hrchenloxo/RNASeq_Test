def calculate_fc(data,change,genes):
    data_ = []
    for no, gene in enumerate(genes):
        data_temp = data[data.index.str.contains(gene)]

        data_temp = data_temp.T

        data_temp['Gene'] = data_temp.index
        #data_temp[['Cell','Day','Condition']] = data_temp.Gene.str.split(r"_", expand=True)[[0,1,2]]
        data_temp[['Cell','Condition','Response']] = data_temp.Gene.str.split(r"_", expand=True)[[0,1,2]]

        ### AUC
        url = "https://raw.githubusercontent.com/hrchenloxo/RNASeq_Test/main/RealAUC_v8.txt"
       # auc = pd.read_table('/Users/hsiao-rongchen/Box/hrchen/Downloads/RXR/Data/RealAUC_v8.txt',header=None)
        auc = pd.read_table(url,header=None)

        auc.columns=['Cell','AUC']
        data_temp = data_temp.merge(auc, left_on='Cell', right_on='Cell',how='left')
        data_temp.sort_values("AUC", inplace=True)
 
        from math import log2

        t = data_temp.groupby(['Cell', 'Condition']).mean()[data_temp.columns[0]]

        comparisons = {}
        for cell, condition in t.index:
            if cell not in comparisons:
                comparisons[cell] = [condition]
            else:
                comparisons[cell].append(condition)

        df_FC = {'Cell': [], 'Compared pair': [], 'FC': []}
        for cell in comparisons:
            if 'Ctrl' not in comparisons[cell]:
                print(f'No control in {cell}. Skip it')
                continue
            ctrl_value = t[(cell, 'Ctrl')]
            for condition in comparisons[cell]:
                if condition == 'Ctrl':
                    continue
                value = t[(cell, condition)]
                if change == 'log2':
                    fc = log2( abs( (value + 0.1)/(ctrl_value + 0.1) ) )
                    #print("Using log2")
                elif change == "diff":
                    fc = (value ) - (ctrl_value)
                    #print("Using difference")
                df_FC['Cell'].append(cell)
                df_FC['Compared pair'].append('%s_Ctrl' % (condition))
                df_FC['FC'].append(fc)
        fc = pd.DataFrame(df_FC)
#         if no == 0:
#             fc.rename(columns={'FC':f'{gene}'},inplace=True)
#         else:
        fc['Gene'] = gene
#         if no == 0:
#             fc_data = fc
      #  else:
            #fc_data = fc_data.merge(fc,left_on=['Cell','Compared pair'],right_on=['Cell','Compared pair'])
        data_.append(fc)

#     if no ==0:
#         fc_data = fc
#     else:
        fc_data = pd.concat(data_)
    return fc_data
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from jupyter_dash import JupyterDash

mergedata1=pd.read_csv("https://raw.githubusercontent.com/hrchenloxo/RNASeq_Test/main/testdata.csv", index_col=0)

fc = calculate_fc(mergedata1, 'log2',[i.split("_")[0] for i in mergedata1.index.unique()][1:5])
#fc = calculate_fc(mergedata1, 'log2',['MYC','RARB'])

JupyterDash.infer_jupyter_proxy_config()
app = JupyterDash(__name__)



df = fc
genes = list(df.Gene.unique())




app.layout = html.Div([
    dcc.Dropdown(
        id="dropdown",
        options=[{"label": x, "value": x} for x in genes],
        value=genes[0],
        clearable=False,
    ),
    dcc.Graph(id="bar-chart"),
])

@app.callback(
    Output("bar-chart", "figure"), 
    [Input("dropdown", "value")])
def update_bar_chart(gene):
    mask = df["Gene"] == gene
    fig = px.bar(df[mask], x="Cell", y="FC", 
                 color="Compared pair", barmode="group")
    return fig

if __name__ == '__main__':
    app.server.run(debug=True, threaded=True)


