# Copyright (c) 2020 Michiel Bongaerts.
# All rights reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
@author: Michiel Bongaerts
"""

import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import seaborn as sns
from collections import defaultdict
import pickle
import networkx as nx
import os
import json


def plot_IEM_ranking_tab(df_ranks, path_to_maps, patient_with_def_gene, block_size = 25):


    all_Ms = [ el.strip('.pickle') for el in os.listdir(path_to_maps)]
    all_Ms_ext = [ el.strip('ext') for el in all_Ms if 'ext' in el]
    all_Ms = [el for el in all_Ms if (el not in all_Ms_ext)]

    def clean_element(x):
        x = x.replace('dr:','_')
        x = x.replace(' rn:','_')
        return x

    path_lengths_R_m = []
    pathway_with_metabs = {}
    for module_or_pathway_ID  in all_Ms:
        if( 'hsa01100' in module_or_pathway_ID):
            continue   

        with open(path_to_maps+'{}.pickle'.format(module_or_pathway_ID.strip('.pickle')), "rb") as f:
            M_edges = pickle.load(f)  

        # Extra stripping on element names
        G_edges = [ ( clean_element(el[0]), clean_element(el[1])  ) for el in M_edges]
        edge_nrm_or_rev = { (el[0],el[1]): el[2] for el in M_edges }

        #make Graph
        G = nx.Graph()
        G.add_edges_from(G_edges)
        G_edges_rev = [(el[1],el[0]) for el in G_edges]
       
        ms = [el for el in G if not (   (  ('R' in el) and ('RM_' not in el )     ) or ('RR_' in el )  ) ]
        Rs = [el for el in G if (   (  ('R' in el) and ('RM_' not in el )     ) or ('RR_' in el )  ) ]
        
        pathway_with_metabs[module_or_pathway_ID] = ms


    df_ranks = df_ranks.assign(relative_rank = df_ranks['rank'] / df_ranks['max_rank'] *100 )

    tab = df_ranks.copy()
    tab['ms'] = tab['ms'].apply(lambda x: json.loads(x.replace("'",'"')) )
    tab['total_measured_metabs'] = tab['ms'].apply(lambda x: len(x))
    tab['total_metabs_in_pathway'] = tab['maps'].apply(lambda x: len(pathway_with_metabs[x]) )

    pat_with_rank = tab.set_index('patient_ID')['rank'].to_dict()
    pat_with_max_rank = tab.set_index('patient_ID')['max_rank'].to_dict()
    pat_with_measured_metabs = tab.set_index('patient_ID')['total_measured_metabs'].to_dict()
    pat_with_metabs_in_pathway = tab.set_index('patient_ID')['total_metabs_in_pathway'].to_dict()


    tab = tab.set_index('patient_ID')[['relative_rank']] #pd.pivot_table( tab, columns = 'comb', index ='patient_ID', values='relative_rank')
    tab = tab.sort_values( by = tab.columns.tolist()[0])
    # tab =tab.assign( _= np.nan )
    tab = tab.assign( **{'Absolute rank': pd.Series(pat_with_rank)})
    tab = tab.assign( **{'Total number of enzymes': pd.Series(pat_with_max_rank)})
    tab = tab.assign( **{'Measured metabolites near reaction': pd.Series(pat_with_measured_metabs)})
    tab = tab.assign( **{'Total metabolites in pathway': pd.Series(pat_with_metabs_in_pathway)})


    sns.set(font_scale=1, style='white')

    tab_renamed = tab.rename(index = {ID: patient_with_def_gene[ID]+' ({})'.format(ID) for ID in tab.index  } )

    cbar_space = 2
    N = (int(tab_renamed.shape[0] / block_size) +1) 

    fig = plt.figure(figsize=( 9*N, 20+ cbar_space ),constrained_layout=True )
    fig.subplots_adjust(wspace=0.001,left=0.15, right=0.95,top=0.975)
    gap = 3
    gs = fig.add_gridspec(block_size +2, N*N + gap*(N-1) )

    for i in range(N):
       
        tab_block = tab_renamed.iloc[(i*block_size): (i+1)*block_size,:]
        tab_block = tab_block.assign( relative_rank = tab_block['relative_rank'].round(2) )
        tab_block = tab_block.assign( **{'Absolute rank' : tab_block['Absolute rank'].round().apply(int)} )
        tab_block = tab_block.assign( **{'Total number of enzymes' : tab_block['Total number of enzymes'].round().apply(int)} )


        if( tab_block.shape[0] < block_size):
            if( tab_block.shape[0]  == 0 ):
                continue
            else:
                ax = fig.add_subplot(gs[cbar_space:tab_block.shape[0], i*N + i*gap:(i+1)*N + (i)*gap ])
                
        else:
            ax = fig.add_subplot(gs[cbar_space:block_size, i*N + i*gap : (i+1)*N + (i)*gap  ])
       
        tab_plot = tab_block.copy()
        # Set color for other cells in table
        tab_plot.iloc[:,1:] = 5
       
        tab_annot = tab_block.applymap(lambda x: str(x) if not np.isnan(x) else '' ).values
        
        tab_plot = tab_plot.rename(columns = {tab_plot.columns[0]: 'Percentile rank'})
               
       
        sns.heatmap(tab_plot, annot = tab_annot, fmt = '', vmax= 10, vmin=0, cmap = 'RdBu_r', ax=ax, 
                    cbar=False,  linewidths=2, linecolor='white', annot_kws = {'fontsize':20} )
                   
       
        ax.set_yticklabels( ax.get_yticklabels(), rotation=0, fontsize=15 )
        ax.set_xticklabels( ax.get_xticklabels(), rotation=90, fontsize=18)
       
        colors  = ['k']
        for xtick, color in zip(ax.get_xticklabels(), colors):
            xtick.set_color(color)


        ax.set_ylabel('')
        ax.set_xlabel('')
       

    cbar_ax = fig.add_subplot(gs[0:cbar_space -1, :   ])

    norm = mpl.colors.Normalize(vmin=0, vmax=10)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap='RdBu_r')
    cmap.set_array([])

    cbar_ax.get_xaxis().set_visible(False)
    cbar_ax.get_yaxis().set_visible(False)
    cbar_ax.axis('off')

    g = fig.colorbar(cmap, ticks= np.arange(0,11,1), ax =cbar_ax, orientation='horizontal', fraction =0.5)

    for t in g.ax.get_xticklabels():
         t.set_fontsize(20)


    cbar_ax.set_title('Percentile rank (%)', fontsize =25)
    plt.show()




def plot_IEM_ranking_curve( df_ranks, name_method = 'Reafect', color = 'red' ):

    df_ranks['PR'] = (df_ranks['rank']/df_ranks['max_rank']).values
    pat_ranks = df_ranks.set_index('patient_ID')['PR'] * 100

    sns.set(font_scale=2.5, style='white')
    fig = plt.figure(figsize=(25,12))
    fig.subplots_adjust(left=0.05, right=0.95,top=0.90, wspace=0.075)

    
    for k,(x_top,y_max,title) in enumerate(zip([101,11],[105,95],['A)','B)']) ):
        ax = fig.add_subplot(1,2,1+k)
        ranks = np.arange(0,x_top -1 + 0.005,0.005)

        top = []
        for rank in ranks:
            top.append(( pat_ranks<= rank).sum())

        y = list(np.array(top)/ len(pat_ranks) *100)
        AUC = int(round(sum(y) * np.diff(ranks)[0]))

        if(x_top < 50 ):
            l = name_method+' partial AUC: '+str(AUC)
        else:
            l = name_method+' AUC: '+str(AUC)

        ax.plot(ranks,y,linewidth=5,
                label= l,
                color = color, linestyle = '-')

        if( k==1):
            for coor in [2.5,5,10]:
                ind_coor = np.where( ranks == coor )[0][0]
                x_coor = ranks[ind_coor]
                y_coor = y[ind_coor]
                ax.scatter([ x_coor ], [ y_coor ], color=color,s=200, edgecolor='k', linewidth=2, zorder=10000)


                ax.text(x_coor + 0.25, y_coor , str(tuple([ round(x_coor,2),int(y_coor)] )),
                        fontsize=30, color = color,
                        bbox= dict(facecolor='w', alpha=0.99))


        ax.set_title(title,fontsize=35, y=1.02)
        if( k == 0):
            ax.set_ylabel('Percentage correct IEM',fontsize=40)
        else:
            ax.set_ylabel('',fontsize=40)
        ax.set_xlabel('Percentile rank',fontsize=35)
        ax.set_ylim([-2,y_max])
        ax.set_xlim([-1,x_top])


        ax.legend(loc=4,fontsize=27, framealpha =1)