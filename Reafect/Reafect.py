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

import os
import sys
import time
import pickle
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import requests
import sympy as sp
from collections import defaultdict
import dask

from .kegg_kgml_parser.KGML_parser import read as read_kgml

## Functions ##

def progressbar(it, prefix="", size=60, file=sys.stdout):
    '''From https://stackoverflow.com/questions/3160699/python-progress-bar '''
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()




## Class ##
class Reafect:
    def __init__(self, path_to_store_maps, path_to_abc_E, path_to_store_results, path_to_databases):

        # Make a path to store all KEGG modules and pathways
        self.path_to_store_maps = path_to_store_maps
        os.makedirs( path_to_store_maps, exist_ok=True)
        
        # Make a path to store all "effective Z-scores". Per module/pathways we determine many paths from metabolite to reaction.
        # These paths are used to calculate the "weighted Z-scores". Since we don't want to redo this calculation all the time, we story all these "wighted Z-scores".
        self.path_to_abc_E = path_to_abc_E
        os.makedirs( path_to_abc_E, exist_ok=True)
        
        # Output directory for results
        self.path_to_store_results = path_to_store_results
        os.makedirs( path_to_store_results, exist_ok=True)
        
        self.path_to_databases = path_to_databases

        # This file needs to be there.
        self.kegg_enzyme_with_reactions = pd.read_csv(path_to_databases+'kegg_enzyme_with_reaction.csv').drop(columns = ['Unnamed: 0'])

        self.decay_paths_abc_last_loaded = (None, None, None)
        
    def get_pathways_kegg(self, hsa_maps):
        '''This function downloads the pathway from http://rest.kegg.jp/ and determines all edges between reaction - and 
        metabolite nodes. Thus, a reaction consists the nodes: substrate --> reaction --> product.
        
        '''
        
        hsa_pathways = {}

        for hsa_map in hsa_maps:
            file = '{}/hsa{}.pickle'.format(  self.path_to_store_maps, hsa_map)
            if( os.path.exists(file)):
                print('{} is already downloaded'.format(hsa_map))
                continue

            # Wait for a request to prevent banishment from server
            time.sleep(0.25)
            print(hsa_map, 'is requested')

            try:
                url  = 'http://rest.kegg.jp/get/hsa{}/kgml'.format(hsa_map)
                res = requests.get(url)
                pathway = read_kgml( res.text )

                edges = [ ]
                for info in pathway.reactions:
                    R = info.name.strip('rn:')+'_'+str(info.id)
                    R = R.replace(' rn:','_')


                    for prod in info.products:
                        edge = (R,prod.name.split(' ')[0].strip('gl:').strip('cpd:').strip('dr:'))

                        if( info.type == 'reversible'):
                            edges.append( (edge[0],edge[1],'nrm') )

                            # Make also the reversed edge, call this 'rv'
                            edges.append( (edge[1],edge[0],'rv') )

                        elif( info.type == 'irreversible'):
                            edges.append( (edge[0],edge[1],'nrm')  )

                        else:
                            print('wrong type')

                    for sub in info.substrates:
                        edge = (sub.name.split(' ')[0].strip('gl:').strip('cpd:').strip('dr:'), R)

                        if( info.type == 'reversible'):
                            edges.append( (edge[0],edge[1],'nrm') )

                            # Make also the reversed edge, call this 'rv'
                            edges.append( (edge[1],edge[0],'rv') )

                        elif( info.type == 'irreversible'):
                            edges.append( (edge[0],edge[1],'nrm') )

                        else:
                            print('Wrong reaction type', info)


                # Remove edges which occur multiple times because they are associated with multiple enzymes 
                # This is necessary to prevent noise in the "weighted Z-scores", else you will get more possible paths.
                df = []
                for edge in edges:
                    a,b,nrm_rv = edge
                    df.append([a,b, (a.split('_')[0], b.split('_')[0], nrm_rv) ])
                df = pd.DataFrame(df, columns = ['a','b','edge_stripped'])


                x = df['edge_stripped'].value_counts()
                df_1 = df.loc[ df['edge_stripped'].isin( x[x==1].index.tolist() )]
                df_2  = df.loc[ df['edge_stripped'].isin(x[x >1].index.tolist())]
                df_2 = df_2.loc[df_2['edge_stripped'].duplicated(keep='first')]
                df = pd.concat([df_1, df_2])

                edges = df['edge_stripped'].tolist()

                # store all edges as a pickle file
                with open( file, 'wb') as handle:
                    pickle.dump( edges, handle, protocol=pickle.HIGHEST_PROTOCOL)

            except Exception as e:
                print(hsa_map, 'Some error has occured:{}'.format(e))

                
                
    def get_info(self, M):
        data_dict = {}
    
        # Wait for a request to prevent banishment from server
        time.sleep(0.25)
        print(M, 'is requested')
        
        result_compound = requests.get('http://rest.kegg.jp/get/{com}'.format(com=M))
        data = [el for el in result_compound.text.split("\n")]
        attribute_list = []
        for row in data:
            attr = row.split(" ")[0]
            if( len(attr.strip("/")) > 1):
                attribute_list.append(attr)

        data_dict = {attr:[] for attr in attribute_list}
        for row in data:
            try:
                attr_ = row.split(" ")[0]
                if attr_ in attribute_list:
                    attr = attr_
                data_dict[attr].append([el for el in row.split("  ")[1:] if el!=''])
            
            except Exception as e:
                print(e)

        return data_dict


    def get_modules_kegg(self, hsa_modules ):
        '''This function downloads the module from http://rest.kegg.jp/ and determines all edges between reaction - and 
        metabolite nodes. Thus, a reaction consists the nodes: substrate --> reaction --> product.
        
        '''
                
        modules = {}
        for M in hsa_modules:

            file = '{}/{}.pickle'.format(self.path_to_store_maps, M)
            if( os.path.exists(file)):
                print('{} is already downloaded'.format(M))
                continue


            module = self.get_info(M)

            if( 'REACTION'  in module):
                edges = []
                Rs = module['REACTION']
                

                for el in Rs:
                    if( len(el) == 0):
                        continue

                    if( len(el) == 1):
                        print("Different input than usual:", el)
                        el_new = []
                        el_new.append(el[0].split(' ')[0])
                        el_new.append(' '.join(el[0].split(' ')[1:]))
                        el = el_new
                        print("Changed input to:",el_new)


                    R = el[0].replace('+','_').replace(',','_')
            

                    for el_ in el:
                        for m in el_.split(','):
                            if( '<->' in m):
                                prods,subs = m.split('<->')
                                if( '+' in prods):
                                    for prod in prods.split("+"):
                                        prod = prod.strip(' ')
                                        edges.append( (prod,R,'nrm' ) )
                                        edges.append( (R ,prod,'rv') )

                                else:
                                    prod = prods.strip(' ')
                                    edges.append( (prod, R, 'nrm' ) )
                                    edges.append( (R, prod, 'rv') )


                                if( '+' in subs):
                                    for sub in subs.split("+"):
                                        sub = sub.strip(' ')
                                        edges.append( (R , sub, 'nrm' ) )
                                        edges.append( (sub, R, 'rv'  ) )

                                else:
                                    sub = subs.strip(' ')
                                    edges.append( (R ,sub, 'nrm' ) )
                                    edges.append( (sub, R, 'rv'  ) )



                            elif( '->' in m):
                                prods,subs = m.split('->')
                                if( '+' in prods):
                                    for prod in prods.split("+"):
                                        prod = prod.strip(' ')
                                        edges.append( (prod, R, 'nrm' ) )

                                else:
                                    prod = prods.strip(' ')
                                    edges.append( (prod, R, 'nrm' ) )

                                if( '+' in subs):
                                    for sub in subs.split("+"):
                                        sub = sub.strip(' ')
                                        edges.append( (R, sub, 'nrm' ) )

                                else:
                                    sub = subs.strip(' ')
                                    edges.append( ( R, sub, 'nrm' ) )

                            else:
                                pass

                # store all edges as a pickle file
                with open(file, 'wb') as handle:
                    pickle.dump( edges, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    
                    
    def update_kegg_pathways(self,extra_reactions, custom_ID_to_KEGG ):
        '''This function opens all stored modules and pathway pickle files located in 'path_to_store_maps' and then adds 
        all the reactions where it finds an overlapping metabolite with a KEGG identifier.
        
        '''
        
        def convert_custom_ID_to_KEGG(x):
            if( x in custom_ID_to_KEGG):
                x = custom_ID_to_KEGG[x]
            return x

        custom_edges = []
        for rev_irrev, path in extra_reactions:
            for i in range(len(path) -1):
                a = convert_custom_ID_to_KEGG( path[i] )
                b = convert_custom_ID_to_KEGG( path[i+1] )

                # Reaction is reversible so we add edges in both directions
                if( rev_irrev == 'rev'):
                    custom_edges.append( (a, b, 'nrm') )
                    custom_edges.append( (b, a, 'rv') )

                elif( rev_irrev == 'irrev'):
                    custom_edges.append( (a, b, 'nrm') )
                else:
                    print('Wrong reaction type',path)

        all_Ms = [ el.strip('.pickle') for el in os.listdir(self.path_to_store_maps)]
        print('amount of modules/pathways', len(all_Ms))
        
        for module_or_pathway_ID  in all_Ms:
            if( 'hsa01100' in module_or_pathway_ID): 
                print(" Map hsa01100 is too big to process ")
                continue            

            # Check if module/pathway already was updated. 
            if( os.path.exists(self.path_to_store_maps+'{}ext.pickle'.format(module_or_pathway_ID) ) ):
                print(" Map {m} has already been updated with additional reactions, delete {m}ext.pickle if you want to redo update.".format(m = module_or_pathway_ID))
                continue
                
            
            # Load edges
            with open(self.path_to_store_maps+'{}.pickle'.format(module_or_pathway_ID), "rb") as f:
                M_edges = pickle.load(f)   
                
            
            # Get all metabolites / reaction nodes from KEGG maps
            nodes = []
            for edge in M_edges:
                a, b, nrm_rv = edge
                nodes.append(a)
                nodes.append(b)

            M_edges = list(set(M_edges))
            N_before = len(M_edges)

            # We loop a few times to make sure all edges are added. Sometimes reactions are added only in the first step
            # but then it turns out another reaction could also 'fit' on the newly added reactions. 
            for i in range(20):

                for edge in custom_edges:
                    a, b, nrm_rv = edge

                    if( a in nodes ):
                        M_edges.append(edge)

                    if( b in nodes ):
                        M_edges.append(edge)

                # update nodes
                nodes = []
                for edge in M_edges:
                    a, b, nrm_rv = edge
                    nodes.append(a)
                    nodes.append(b)

            M_edges = list(set(M_edges))
            N_after = len(M_edges)
            
            if( (N_before - N_after) != 0 ):
                print( "Module/pathway {} has been updated. #edges before {}, #edges after {}.".format(module_or_pathway_ID, N_before, N_after ) )

                with open(self.path_to_store_maps+'{}ext.pickle'.format( module_or_pathway_ID ), 'wb') as handle:
                    pickle.dump( M_edges, handle, protocol=pickle.HIGHEST_PROTOCOL)


                                                

    def prepare_pathways(self):

        file = '{}/all_m_R_E.pickle'.format(self.path_to_abc_E)
        if( os.path.exists(file)):
            print('Pathways already prepared. Delete {} if you want to redo pathway preparation. You should do this when you manually added more reactions.'.format(file))
            
        else:

            def clean_element(x):
                x = x.replace('dr:','')
                x = x.replace(' rn:','_')
                return x
            
            # Make sympy variables for decay factors 
            all_vars = ['a','b','c']
            sym_vars = {}
            for var in all_vars:
                sym_vars[var] = sp.symbols(var)

                              #Z-sign  #upstream           downstream         # reversible reaction
            loss_fractions = { 1: {1 : sym_vars['a'], -1 : sym_vars['b'], 0:sym_vars['c']}, # pos Z-score
                              -1: {1 : sym_vars['b'], -1 : sym_vars['a'], 0:sym_vars['c']},  # neg Z-score
                             }    

            # Get all pathways and modules 
            all_Ms = [ el.strip('.pickle') for el in os.listdir(self.path_to_store_maps)]
            # Remove all original maps/pathways which are extended with reactions
            all_Ms_ext = [ el.strip('ext') for el in all_Ms if 'ext' in el]
            all_Ms = [el for el in all_Ms if (el not in all_Ms_ext)]


            decay_paths = []
            for module_or_pathway_ID  in all_Ms:
                print(module_or_pathway_ID)

                if( 'hsa01100' in module_or_pathway_ID): 
                    print(" Map hsa01100 is too big to process ")
                    continue    

                with open('{}/{}.pickle'.format(self.path_to_store_maps, module_or_pathway_ID.strip('.pickle')), "rb") as f:
                    M_edges = pickle.load(f)   


                # Extra stripping on element names
                G_edges = [ ( clean_element(el[0]), clean_element(el[1])  ) for el in M_edges]
                edge_nrm_or_rev = { (el[0],el[1]): el[2] for el in M_edges }

                # Make Graph
                G = nx.Graph()
                G.add_edges_from(G_edges)

                # Get all metabolite nodes in the graph
                ms = [el for el in G if not (   (  ('R' in el) and ('RM_' not in el )     ) or ('RR_' in el )  ) ]
                # Get all reaction nodes in the graph
                Rs = [el for el in G if (   (  ('R' in el) and ('RM_' not in el )     ) or ('RR_' in el )  ) ] 

                # Make sympy variables for each metabolite
                for m in ms:
                    if( m not in sym_vars):
                        sym_vars[m] = sp.symbols(m)

                # Iterate over all metabolites
                for m in ms:
                    # Select network in the neighbourhood of metabolite m
                    ego = nx.ego_graph(G, m, radius = 15)
                    Rs_ego = [el for el in ego if (   (  ('R' in el) and ('RM_' not in el )     ) or ('RR_' in el )  ) ] 
                    # Determine all paths for m to all reactions
                    all_pairs = list(itertools.product([m], Rs_ego))


                    # Iterate over all metabolite-reaction pairs in ego-graph
                    for pair in all_pairs:
                        (s,t) = pair
                        path_counter = 0
                        
                        all_paths = []
                        
                        # Iterate over all cuttoffs, make sure we start with lower cutoffs because these are most important
                        # i.e. most near metabolites
                        for cuttoff in [3,5,7,9,11,13,15]:
                            # Search all simple paths from metabolite to reaction
                            for ii,path in enumerate(nx.all_simple_paths(ego, source  = s, target  = t, cutoff = cuttoff)):
                                
                                # Don't continue if path already has been used
                                if( tuple(path) in all_paths):
                                    continue
                                
                                # Don't continue if max number of paths has been reached
                                if( path_counter >  10 ):
                                    break

                                edges = [(path[i],path[i+1]) for i in range(len(path)-1)]

                                sign_stream = []
                                rev_irrev = []
                                decay_pos  = []
                                decay_neg  = []
                                for edge in edges:

                                    # # Is the edge from a reversible reaction ? Then both directions will exists in G_edges
                                    if( (edge in G_edges) and ((edge[1],edge[0]) in G_edges) ):

                                        decay_pos.append( loss_fractions[1][0] )
                                        decay_neg.append( loss_fractions[-1][0] )

                                        # We take the definition of upstream or downstream for reversible reaction by
                                        # how they were added. 'nrm' edges were added first and 'rv' edges where the 'nrm' edges reversed.
                                        # So we take 'nrm' as the definition of up- / downstream. 
                                        if( edge_nrm_or_rev[edge] == 'nrm' ):
                                            sign_stream.append(1)
                                            rev_irrev.append('rev')
                                        else:
                                            sign_stream.append(-1)
                                            rev_irrev.append('rev')

                                    # Downstream
                                    elif(edge in G_edges):
                                        decay_pos.append( loss_fractions[1][1] )
                                        decay_neg.append( loss_fractions[-1][1] )
                                        sign_stream.append(1)
                                        rev_irrev.append('irrev')

                                    # upstream
                                    elif( (edge[1],edge[0]) in G_edges):
                                        decay_pos.append( loss_fractions[1][-1] )
                                        decay_neg.append( loss_fractions[-1][-1] )
                                        sign_stream.append(-1)
                                        rev_irrev.append('irrev')

                                    # unknown
                                    else:
                                        pass

                                path_counter+=1
                                
                                # Append path 
                                all_paths.append( tuple(path) )

                                decay_paths.append(  { 'decay_pos':np.prod(decay_pos) * sym_vars[s],
                                                       'decay_neg':np.prod(decay_neg) * sym_vars[s],
                                                       'last_step_direct':sign_stream[-1], # This indicates if the connection is on the up- or downstream side of the reaction
                                                       'rev_irrev':rev_irrev[-1], 
                                                       'm': s, # metabolite
                                                       'R': t, # reaction node name
                                                       'maps': module_or_pathway_ID, #module/pathway name
                                                       'p': tuple(path),
                                                      })
                                
            # Store these results as pickle. We do not want to calculate this twice.
            with open(file, 'wb') as handle:
                pickle.dump(decay_paths, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
            
            
            
    def determine_pathways_for_abc_params(self, a_params = [0.85] , b_params = [0.35], c_params = [0.85], n_cores = 4 ):

        # Check input of a,b,c
        cnt =  sum([ (a>1) or (a<0) for a in a_params])
        cnt += sum([ (b>1) or (b<0) for b in b_params])
        cnt += sum([ (c>1) or (c<0) for c in c_params])

        if( cnt > 0):
            raise ValueError("Wrong input for a_params, b_params, c_params. Numbers need to be between 0 and 1.")

        
        with open(self.path_to_abc_E+'all_m_R_E.pickle', "rb") as f:
            decay_paths = pickle.load(f)   
        df_decay_paths = pd.DataFrame(decay_paths)


        def subs_abc(a_param,b_param,c_param):
            print( 'a={}, b={}, c={}'.format(a_param,b_param,c_param) )

            def subs_ab(x):
                try:
                    return x.subs({'a':a_param,'b':b_param,'c':c_param})
                except:
                    return x

            df_decay_paths_abc = df_decay_paths.copy()
            df_decay_paths_abc['decay_pos'] = df_decay_paths_abc['decay_pos'].apply(subs_ab)
            df_decay_paths_abc['decay_neg'] = df_decay_paths_abc['decay_neg'].apply(subs_ab)

            with open(self.path_to_abc_E+'all_m_R_E_{a}_{b}_{c}.pickle'.format(a=a_param,b=b_param,c=c_param), 'wb') as handle:
                pickle.dump(df_decay_paths_abc, handle, protocol=pickle.HIGHEST_PROTOCOL)


        outputs = [] 
        # Go over all given numbers of a,b,c
        for a_param in a_params:
            for b_param in b_params:
                for c_param in c_params:
                    file  = self.path_to_abc_E+'all_m_R_E_{a}_{b}_{c}.pickle'.format(a=a_param,b=b_param,c=c_param)
                    
                    # Check if pathway for this a,b,c has been made already.
                    if( not os.path.exists(file) ):
                        outputs.append( dask.delayed(subs_abc)( a_param,b_param,c_param )  )

        # We just make a pseudo function used for parallel computing
        def concat(outputs):
            return 1

        _ = dask.delayed(concat)(outputs).compute( num_workers = n_cores)


        
    def process_patient(self, patient_ID,  Z_score_data, type_of_input_Z_scores = 'Z_scores', disease_enzyme = None,
                        a_param = 0.85 , b_param = 0.35, c_param = 0.75, custom_reactions_with_enzyme = pd.DataFrame([], columns = ['enzyme','R'] )  ):

        # Check input of custom_reactions_with_enzyme
        assert isinstance(custom_reactions_with_enzyme, pd.DataFrame )
        assert( len(set(custom_reactions_with_enzyme.columns).intersection(['enzyme','R'])) == 2 )

        def subs_Z_score(row):
            # Z-score is positive, so all decay paths for negative Z-scores are set to nan
            if( np.sign(row['Z_score']) == 1):
                row['decay_pos'] = row['decay_pos'].subs({ row['m']: row['Z_score'] })
                row['decay_neg'] = np.nan

            # Z-score is negative, so all decay paths for positive Z-scores are set to nan
            elif( np.sign(row['Z_score']) == -1):
                row['decay_neg'] = row['decay_neg'].subs({ row['m']: row['Z_score'] })
                row['decay_pos'] = np.nan

            else:
                row['decay_neg'] = np.nan
                row['decay_pos'] = np.nan

            return row


        print("Started", patient_ID )
        path_to_sample = self.path_to_store_results+'{}/'.format(patient_ID)
        os.makedirs( path_to_sample, exist_ok=True)


        file = path_to_sample+'SR_scores_{a}_{b}_{c}_{type_of_input_Z_scores}.csv'.format(a=a_param, b=b_param, c=c_param, type_of_input_Z_scores = type_of_input_Z_scores)
        if( not os.path.exists(file) ):

            if( self.decay_paths_abc_last_loaded != (a_param, b_param, c_param) ):
                with open(self.path_to_abc_E+'all_m_R_E_{a}_{b}_{c}.pickle'.format(a = a_param, b = b_param, c = c_param ), "rb") as f:
                    self.df_decay_paths_abc = pickle.load(f)   
                    self.decay_paths_abc_last_loaded = (a_param, b_param, c_param)

            df_pat_abc = self.df_decay_paths_abc.merge( Z_score_data.reset_index().rename(columns = {'network_ID': 'm'}), on ='m' )


            SR_scores = []
            gb = df_pat_abc.groupby(by=['maps','R'])
            for i,( ( ( module_pathway, R ), data_R), progress) in enumerate(zip(gb, progressbar(range(len(gb)), "Computing: ", 60) ) ):


                data_R = data_R.apply( subs_Z_score, axis=1)

                data_R_pos = data_R.loc[ data_R['decay_pos'].notnull() ]
                if( not data_R_pos.empty):
                    # Group by metabolite and calculate normalized Z-score for every path
                    weighted_E_pos = data_R_pos.groupby(by = ['m'] ).apply(lambda x: x[['decay_pos']]**2 / x[['decay_pos']].abs().sum()   )
                    weighted_E_pos = weighted_E_pos.rename(columns = {'decay_pos' : 'decay_pos_weighted'})
                    data_R = data_R.assign(decay_pos_weighted = weighted_E_pos['decay_pos_weighted'])
                else:
                    data_R['decay_pos_weighted'] = np.nan

                data_R_neg = data_R.loc[ data_R['decay_neg'].notnull() ]
                if( not data_R_neg.empty):
                    # Group by metabolite and calculate normalized Z-score for every path
                    weighted_E_neg = data_R_neg.groupby(by = ['m'] ).apply(lambda x: x[['decay_neg']]**2 / x[['decay_neg']].abs().sum()   )
                    weighted_E_neg = weighted_E_neg.rename(columns = {'decay_neg' : 'decay_neg_weighted'})
                    data_R = data_R.assign(decay_neg_weighted = weighted_E_neg['decay_neg_weighted'])

                else:
                    data_R['decay_neg_weighted'] = np.nan

                SR = 0
                for up_down, data_up_down in data_R.groupby('last_step_direct'):
                    # upstream side of reaction
                    if( up_down == 1):
                        SR_up = sum(data_up_down['decay_pos_weighted'].dropna()) - sum(data_up_down['decay_neg_weighted'].dropna())
                        SR+=SR_up

                    # downstream side of reaction
                    if( up_down == -1):
                        SR_down = sum(data_up_down['decay_neg_weighted'].dropna()) - sum(data_up_down['decay_pos_weighted'].dropna())  
                        SR+=SR_down

                
                uni_rev_irrev = data_R['rev_irrev'].unique().tolist()

                # Aparrently this reaction is neither reversible or irreversible
                # Then you should check what went wrong!
                if( len(uni_rev_irrev) > 1):
                    print( "Reaction {} has no clear directionality!".format(R) )
                    SR_scores.append([R,SR, list(data_R['m'].unique() ), module_pathway, uni_rev_irrev ])

                elif( uni_rev_irrev[0] == 'irrev'):
                    SR_scores.append([R, SR, list(data_R['m'].unique() ), module_pathway,'irrev' ])

                # reversible reaction than absolute value of SR
                elif( uni_rev_irrev[0] == 'rev'):
                    SR_scores.append([R, abs(SR), list(data_R['m'].unique() ), module_pathway, 'rev' ])

                # Check what went wrong!
                else:
                    SR_scores.append([R, SR, list(data_R['m'].unique() ), module_pathway, 'error' ])


            SR_scores  = pd.DataFrame(SR_scores, columns = ['R','SR','ms','maps','irrev_rev'])
            SR_scores['SR'] = SR_scores['SR'].apply(np.float64)
            SR_scores = SR_scores.sort_values(by='SR', ascending=False).reset_index(drop=True)
            SR_scores['disease_enzyme'] = disease_enzyme
            SR_scores['patient_ID'] = patient_ID
            
            # Copy enzyme with reaction database but also add user input additional associations
            kegg_enzyme_with_reactions = self.kegg_enzyme_with_reactions.copy()
            kegg_enzyme_with_reactions = pd.concat([kegg_enzyme_with_reactions, custom_reactions_with_enzyme])
            kegg_enzyme_with_reactions = kegg_enzyme_with_reactions.loc[ ~ kegg_enzyme_with_reactions.duplicated(keep='first') ]
            kegg_enzyme_with_reactions = kegg_enzyme_with_reactions.groupby('enzyme').apply(lambda x: x['R'].tolist()).to_dict()


            # Search enzymes with reactions
            enzyme_with_R = []
            unique_R = SR_scores.loc[ ~SR_scores['R'].duplicated(keep='first')]

            for enz,Rs in kegg_enzyme_with_reactions.items():
                
                # make string with all reaction belonging to an enzyme
                s = ''
                for R in Rs:
                    s += R+'|'
                s = s.strip('|')

                if( len(s) > 4):
                    # Search matches 
                    R_matches = unique_R.loc[ unique_R['R'].str.contains(s),'R'].tolist()
                    for R in R_matches:
                        enzyme_with_R.append([enz,R])
                        
            
            enzyme_with_R = pd.DataFrame(enzyme_with_R, columns = ['enzyme','R'])
            enzyme_with_R = enzyme_with_R.loc[ ~ enzyme_with_R.duplicated(keep='first') ]
            SR_scores_enz = SR_scores.merge(enzyme_with_R, on='R', how = 'left')

            SR_scores_enz = SR_scores_enz.sort_values(by= ['SR','enzyme'],ascending=False).reset_index(drop=True)
            ranks = {}
            rank = 0
            for enz in SR_scores_enz['enzyme'].dropna().tolist():
                if( enz in ranks):
                    pass
                else:
                    ranks[enz] = rank
                    rank+=1

            SR_scores_enz['rank'] = SR_scores_enz['enzyme'].map(ranks)
            SR_scores_enz['max_rank'] = SR_scores_enz['rank'].max()

            SR_scores_enz.to_csv(file)
            print("Done and results are stored at {}".format(file) )

            return SR_scores_enz

        else:
            SR_scores_enz = pd.read_csv(file)
            print("Results are already determined and stored. Delete file if you want to redo processing.")
            return SR_scores_enz

        
        

        