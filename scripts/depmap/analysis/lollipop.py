import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

cwd = os.getcwd()  

file_name = "analysis/test_lollipop/RAD51C_clusters_methylation_STAD.csv"
csv_location = os.path.join(cwd, file_name)

df = pd.read_csv(csv_location)
# Rename columns without overwriting original position values
column_names = ["cell_line"] + [str(i) for i in df.columns[1:-1]] + ["level"]
df.columns = column_names

print(column_names)

# Extract the specified cell lines
#selected_cell_lines = ["SNU601", "AGS", "TGBC11TKB", "SNU620", "SNU719", "SH10TC", "GSS", "LMSU", "MKN74"]
#selected_cell_lines = ["SNU601", "SNU620", "SNU719", "SH10TC", "GSS"]
#WXS cell lines only
#selected_cell_lines = ['GSS','GSU','HUG1N','LMSU','MKN74','NCCSTCK140','NUGC2','NUGC4','SH10TC','SNU216','SNU520','SNU601','SNU620','SNU668','SNU719']
selected_cell_lines = df.cell_line.to_list()
print(selected_cell_lines)

df_selected = df[df["cell_line"].isin(selected_cell_lines)]


# Prepare the dataframe for plotting
df_melt = df_selected.melt(id_vars="cell_line", value_vars=column_names[1:-1])
df_melt.columns = ['variable', 'pos', 'value']
# using the RAD51C exon start site 
df_melt['pos'] = df_melt['pos'].astype(int) - 56770005
#add two to position if it equals 76 to space out points
df_melt.loc[df_melt['pos'] == 76, 'pos'] += 3
# set uniform space between cell lines across panels
cell_line_positions = {cell_line: idx for idx, cell_line in enumerate(selected_cell_lines)}

#classifications = ["no methylation", "homozygous", "heterozygous"]
classifications = ["homozygous", "heterozygous", "no methylation"]

yticks = {}
for classification in classifications:
    yticks[classification] = df_selected[df_selected["level"] == classification]["cell_line"].tolist()

max_samples_in_group = max(df_selected['level'].value_counts())  

# adjust subplot heights based on the number of samples
#height_ratios = [len(yticks[classification]) for classification in classifications]
height_ratios = [0.05, 0.2, 1.35]
print(height_ratios)

title_mapping = {
    "homozygous": "Homozygous meRAD51C",
    "heterozygous": "Heterozygous meRAD51C",
    "no methylation": "No Methylation"
}

fig, ax = plt.subplots(3, 1, figsize=(8, 11), sharex=True, gridspec_kw={'height_ratios': height_ratios})

for idx, classification in enumerate(classifications):
    subset = df_melt[df_melt["variable"].isin(yticks[classification])]
    
    Y_TICK_FRACTION = 0.65
    # mapping for y-values within each subplot
    #local_y_mapping = {cell_line: i for i, cell_line in enumerate(subset['variable'].unique())}
    local_y_mapping = {cell_line: (i * Y_TICK_FRACTION) + 0.05 for i, cell_line in enumerate(subset['variable'].unique())}


    subset = subset.copy()
    subset['local_y'] = subset['variable'].map(local_y_mapping)
    
    ax[idx].hlines(subset['local_y'], min(subset['pos']), max(subset['pos']), label='_nolegend_', zorder=1)
    
    # remove edgecolors for scatter points
    im = ax[idx].scatter(subset['pos'], subset['local_y'], c=subset['value'], s=100,
                         edgecolors="black", cmap='coolwarm', zorder=2)
    im.set_clim(0, 1)

    # remove tick marks
    ax[idx].tick_params(axis='both', which='both', length=0)
    ax[idx].set_yticks(list(local_y_mapping.values()))
    
    #ax[idx].set_ylim(-1, len(local_y_mapping))
    if len(local_y_mapping) > 1:
        specific_y_lim = len(local_y_mapping)* Y_TICK_FRACTION 
    else:
        specific_y_lim = 0.3

    ax[idx].set_ylim(-0.25, specific_y_lim)  
    ax[idx].set_yticklabels(list(local_y_mapping.keys()))

    if idx == 0: 
        specific_y = 0.7
    elif idx == 1:
        specific_y = 0.8
    else: 
        #could be made proportional to the number of samples for this plot
        specific_y = 0.97

    #specific_y = height_ratios[idx]+ 1/(idx + 1)


    # set the subplot title (label) to the side
    ax[idx].set_title(title_mapping[classification], loc='right', fontsize= 10, fontweight='bold',  y=specific_y)

    # remove borders from the subplot
    for spine in ax[idx].spines.values():
        spine.set_visible(False)

fig.subplots_adjust(hspace=0.05)  # adjust 0.2 as needed

# colorbar and x-axis 
cbar = fig.colorbar(im, ax=ax.ravel().tolist(), ticks=[0.1, 0.5, 0.9], orientation='horizontal', 
                    pad=0.05,  # Increase padding for moving it further down 
                    shrink=0.6,  # Adjust the shrink value as needed to change the size
                    aspect=30  # Adjust the aspect ratio as needed
                   )
cbar.ax.tick_params(labelsize=6.5, length=1.5, pad=0.2)
cbar.set_label("CpG methylation", size=7.5)  # Add title for the colorbar


ax[-1].set_xticks(df_melt['pos'].unique())
ax[-1].tick_params(axis='x', which='major', labelsize=7, rotation=45, pad = 0.2)
ax[-1].set_xlabel("meRAD51C CpG site", fontsize=8, labelpad = 0.05)  # Add x-axis title

# Change the tick label
xtick_labels = ax[-1].get_xticks().tolist()
xtick_labels = [str(int(label)) if label != 79 else '76' for label in xtick_labels]
ax[-1].set_xticklabels(xtick_labels)

outpath = os.path.join(cwd, "analysis/test_lollipop/")

outname="All_samples_combined_colour_facet.pdf"

#plt.tight_layout()
fig.savefig(outpath + "/" + outname)
plt.close()
#plt.show()