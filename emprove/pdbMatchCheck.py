#!/usr/bin/python3.8

# !/usr/local/bin/python3.9


from math import sqrt
from Bio.PDB import PDBParser, PDBIO, Chain, Residue, Model, Structure #, PDBConstructionException
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os,sys
import seaborn as sns
import mrcfile

__version__ = '0.13'



#### COMMON FUNCTIONS ################
def calc_distance(coord1, coord2):
    """
    Calculate the Euclidean distance between two points in 3D.
    Returns a float representing the distance.
    """
    distance = sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)
    return distance

def count_matches(df_group, name_column1, name_column2):
    """
    Count the occurrences where residue1 is equal to residue2 in a DataFrame group.
    Returns a dictionary with the count.
    """
    match_count = (df_group[name_column1] == df_group[name_column2]).sum()
    return {'match': match_count}
##########################################

#compute Residues List
def getResiduesList(filename, mask_file=None):
    """
    Compute and return a list of residues from a PDB structure file.
    Each residue in the list is a tuple including residue name, average coordinates of its atoms, chain ID, and residue ID.
    Only residues corresponding to the mask are included in the list if a mask file is provided.
    """
    if mask_file is not None:
        with mrcfile.open(mask_file, 'r') as mrc:
            numpy_image = mrc.data.copy() # make a copy to avoid modifying the original data
            voxel_size = np.array([mrc.voxel_size.x.item(), mrc.voxel_size.y.item(), mrc.voxel_size.z.item()])
            origin = np.array([mrc.header['origin']['x'].item(), mrc.header['origin']['y'].item(), mrc.header['origin']['z'].item()])
        numpy_image = np.transpose(numpy_image, (2, 1, 0))
    else:
        numpy_image = None


    print('residues reading')
    parser = PDBParser()
    try:
        model1 = parser.get_structure('model1', filename)[0]
    except:
        print('    some errors happens reading ', filename)

    residue_list = []
    # Iterate over all the chains in the structure
    for chain in model1.get_chains():
        for residue in chain:
            residue_name = residue.get_resname()
            residue_atoms = residue.get_atoms()
            residue_id = residue.get_id()
            num_residues = 0
            coords= np.array([0.0, 0.0, 0.0])
            #residueIDs=[]
            for atom in residue_atoms:
                #if mask_file is not None:
                #    voxel_idx = np.round((atom.coord - origin) / voxel_size).astype(int)
                #    if np.all(0 <= voxel_idx) and np.all(voxel_idx < numpy_image.shape):
                #        if numpy_image[tuple(voxel_idx)] <= 0.1:
                #            continue
                num_residues+=1
                coords+=atom.get_coord()
                #residueIDs.append(atom.get_serial_number())
            if num_residues > 0:  # Only include residues that have at least one atom in the mask
                coords=coords/num_residues
                #float_list = [float(elem) for elem in coords]
                residue_list.append((residue_name, coords, chain.get_id() ,residue.get_id()[1]))
                #print (residue_name,"    chainID=", chain.get_id(),"   ID=",residue.get_id()[1])
    return residue_list



#compute Full Atoms List (Not used)
def getAtomsList(filename, mask_file=None):
    """
    Compute and return a list of all atoms from a PDB structure file.
    Each atom in the list is a tuple including atom name, its coordinates, chain ID, and residue ID.
    Only atoms corresponding to the mask are included in the list if a mask file is provided.
    """
    if mask_file is not None:
        with mrcfile.open(mask_file, 'r') as mrc:
            numpy_image = mrc.data.copy() # make a copy to avoid modifying the original data
            voxel_size = np.array([mrc.voxel_size.x.item(), mrc.voxel_size.y.item(), mrc.voxel_size.z.item()])
            origin = np.array([mrc.header['origin']['x'].item(), mrc.header['origin']['y'].item(), mrc.header['origin']['z'].item()])
        numpy_image = np.transpose(numpy_image, (2, 1, 0))
    else:
        numpy_image = None


    parser = PDBParser()
    try:
        model1 = parser.get_structure('model1', filename)[0]
    except:
        print('    some errors happens opening ',filename)

    atom_list = []
    # Iterate over all the chains in the structure
    for chain in model1.get_chains():
        for residue in chain:
            residue_atoms = residue.get_atoms()
            for atom in residue_atoms:
                print ('name=',atom.get_name())
                #        if atomList1[ii][0]=="CA" or atomList1[ii][0]=="N" or atomList1[ii][0]=="C":
                if mask_file is not None:
                    voxel_idx = np.round((atom.coord - origin) / voxel_size).astype(int)
                    if np.all(0 <= voxel_idx) and np.all(voxel_idx < numpy_image.shape):
                        if numpy_image[tuple(voxel_idx)] <= 0.1:
                            continue
                coords=atom.get_coord()
                atom_list.append((atom.get_name(), coords, chain.get_id() ,residue.get_id()[1]))
    return atom_list

#compute Backbone List
def getBackboneAtomsList(filename, mask_file=None):
    """
    Compute and return a list of backbone atoms (N, C, CA) from a PDB structure file.
    Each atom in the list is a tuple including atom name, its coordinates, chain ID, and residue ID.
    Only atoms corresponding to the mask are included in the list if a mask file is provided.
    """
    if mask_file is not None:
        with mrcfile.open(mask_file, 'r') as mrc:
            numpy_image = mrc.data.copy() # make a copy to avoid modifying the original data
            voxel_size = np.array([mrc.voxel_size.x.item(), mrc.voxel_size.y.item(), mrc.voxel_size.z.item()])
            origin = np.array([mrc.header['origin']['x'].item(), mrc.header['origin']['y'].item(), mrc.header['origin']['z'].item()])
        numpy_image = np.transpose(numpy_image, (2, 1, 0))
    else:
        numpy_image = None

    print('backbone reading')
    parser = PDBParser()
    try:
        model1 = parser.get_structure('model1', filename)[0]
    except:
        print('    some errors happens reading', filename)

    atom_list = []
    # Iterate over all the chains in the structure
    for chain in model1.get_chains():
        for residue in chain:
            residue_atoms = residue.get_atoms()
            for atom in residue_atoms:
                if atom.get_name() == "CA" or atom.get_name()=="N" or atom.get_name()=="C":
                    if mask_file is not None:
                        voxel_idx = np.round((atom.coord - origin) / voxel_size).astype(int)
                        if np.all(0 <= voxel_idx) and np.all(voxel_idx < numpy_image.shape):
                            if numpy_image[tuple(voxel_idx)] <= 0.1:
                                continue
                    coords=atom.get_coord()
                    atom_list.append((atom.get_name(), coords, chain.get_id() ,residue.get_id()[1]))
    return atom_list


def maskPdbSelect(pdbFileIn, pdbFileOut, mask_file):
    import mrcfile
    with mrcfile.open(mask_file, 'r') as mrc:
        numpy_image = mrc.data.copy() # make a copy to avoid modifying the original data
        voxel_size = np.array([mrc.voxel_size.x.item(), mrc.voxel_size.y.item(), mrc.voxel_size.z.item()])
        origin = np.array([mrc.header['origin']['x'].item(), mrc.header['origin']['y'].item(), mrc.header['origin']['z'].item()])
    numpy_image = np.transpose(numpy_image, (2, 1, 0))

    print ('numpy_image.shape=',numpy_image.shape)
    parser = PDBParser()
    try:
        structure = parser.get_structure('model', pdbFileIn)[0]
    except:
        print('    some errors happens opening ', pdbFileIn)

    # Prepare to write a new PDB file
    io = PDBIO()

    # Filter residues and write output PDB
    accepted_residues = []
    for residue in structure.get_residues():
        for atom in residue:
            voxel_idx = np.round((atom.coord - origin) / voxel_size).astype(int)
            if np.all(0 <= voxel_idx) and np.all(voxel_idx < numpy_image.shape):
                if numpy_image[tuple(voxel_idx)] > 0.1:
                    accepted_residues.append(residue)
                    break

    # New structure only containing accepted residues
    new_structure = Structure.Structure('new_structure')
    new_model = Model.Model(0)
    new_chain = Chain.Chain('A')
    for i, residue in enumerate(accepted_residues, start=1):
        new_residue = Residue.Residue(residue.get_id(), residue.get_resname(), residue.get_segid())
        for atom in residue:
            new_residue.add(atom)
        new_chain.add(new_residue)
    new_model.add(new_chain)
    new_structure.add(new_model)

    io.set_structure(new_structure)
    io.save(pdbFileOut)



def getResiduesCorrespondenceMatrix(residuesList1, residuesList2):
    """
    Compute and return a DataFrame containing a correspondence matrix between two sets of residues based on the spatial proximity.
    Each row of the DataFrame includes information about a pair of closest residues from two structures.
    """
    df = pd.DataFrame(columns=['residue1', 'residue2', 'chainResidue1', 'idResidue1', 'chainResidue2', 'idResidue2', 'distance'])

    for ii in range(0,len(residuesList1)):
        target=0
        distance = 9999999        
        for jj in range(0,len(residuesList2)):
            distanceTmp=calc_distance(residuesList1[ii][1], residuesList2[jj][1])
            if distanceTmp < distance:
                distance=distanceTmp
                target=jj

        new_row = {'residue1': residuesList1[ii][0], 'residue2': residuesList2[target][0],   'chainResidue1': residuesList1[ii][2], 'idResidue1': residuesList1[ii][3],  'chainResidue2': residuesList2[target][2], 'idResidue2': residuesList2[target][3], 'distance': distance}
        df = pd.concat([df, pd.DataFrame(new_row, index=[0])], ignore_index=True)
        
        #print (new_row)

        #df = df.append({'residue1': residuesList1[ii][0], 'residue2': residuesList2[target][0], 'chainResidue1': residuesList1[ii][2], 'idResidue1': residuesList1[ii][3], 'chainResidue2': residuesList2[target][2], 'idResidue2': residuesList2[target][3], 'distance': distance}, ignore_index=True)
    return df



def getBackboneCorrespondenceMatrix(atomList1, atomListGT):
    """
    Compute and return a DataFrame containing a correspondence matrix between two sets of backbone atoms based on the spatial proximity.
    Each row of the DataFrame includes information about a pair of closest backbone atoms from two structures.
    """
    df = pd.DataFrame(columns=['atom1', 'atom2', 'chainResidue1', 'idResidue1', 'chainResidue2', 'idResidue2', 'distance'])
    for ii in range(0,len(atomList1)):
        if atomList1[ii][0]=="CA" or atomList1[ii][0]=="N" or atomList1[ii][0]=="C":
            #print (atomList1[ii])
            target=0
            distance = 9999999
            print 
            for jj in range(0,len(atomListGT)):
                distanceTmp=calc_distance(atomList1[ii][1], atomListGT[jj][1])
                if distanceTmp < distance:
                    distance=distanceTmp
                    target=jj

            new_row = {'atom1': atomList1[ii][0], 'atom2': atomListGT[target][0], 'chainResidue1': atomList1[ii][2], 'idResidue1': atomList1[ii][3], 'chainResidue2': atomListGT[target][2], 'idResidue2': atomListGT[target][3], 'distance': distance}
            df = pd.concat([df, pd.DataFrame(new_row, index=[0])], ignore_index=True)

            #df = df.append({'atom1': atomList1[ii][0], 'atom2': atomList2[target][0], 'chainResidue1': atomList1[ii][2], 'idResidue1': atomList1[ii][3], 'chainResidue2': atomList2[target][2], 'idResidue2': atomList2[target][3], 'distance': distance}, ignore_index=True)
    return df


def identify_mismatches_and_missing(df1, df2, thresholdDistance):
    """
    Identify exact matches, incorrect or missing residues from both DataFrames.
    """
    exact_matches = []
    incorrect_or_missing_df1 = []
    incorrect_or_missing_df2 = []
    
    # Identify exact matches and mismatches in df1
    for index, row in df1.iterrows():
        if row['distance'] <= thresholdDistance:
            if row['residue1'] == row['residue2']:
                exact_matches.append((row['residue1'], row['chainResidue1'], row['idResidue1']))
            else:
                incorrect_or_missing_df1.append((row['residue1'], row['chainResidue1'], row['idResidue1']))
        else:
            incorrect_or_missing_df1.append((row['residue1'], row['chainResidue1'], row['idResidue1']))
    
    # Identify mismatches in df2 (incorrect or missing in df1)
    for index, row in df2.iterrows():
        if row['distance'] <= thresholdDistance and row['residue1'] != row['residue2']:
            incorrect_or_missing_df2.append((row['residue2'], row['chainResidue2'], row['idResidue2']))
        elif row['distance'] > thresholdDistance:
            incorrect_or_missing_df2.append((row['residue2'], row['chainResidue2'], row['idResidue2']))

    # Convert to DataFrame for easier handling
    results = {
        'Exact Matches': exact_matches,
        'Incorrect/Missing from DF1': incorrect_or_missing_df1,
        'Incorrect/Missing from DF2': incorrect_or_missing_df2,
    }
    return results


def plot_hits_and_misses(results):
    """
    Plot aligned bar charts for hits and misses based on the comparison results.
    
    Parameters:
    - results: A dictionary containing the comparison results from identify_mismatches_and_missing function.
    """
    # Extracting residue names and their counts for hits, misses from DF1, and misses from DF2
    exact_matches_names = [res[0] for res in results['Exact Matches']]
    missing_from_df1_names = [res[0] for res in results['Incorrect/Missing from DF1']]
    missing_from_df2_names = [res[0] for res in results['Incorrect/Missing from DF2']]
    
    # Counting occurrences
    hits_counts = {name: exact_matches_names.count(name) for name in set(exact_matches_names)}
    missing_df1_counts = {name: missing_from_df1_names.count(name) for name in set(missing_from_df1_names)}
    missing_df2_counts = {name: missing_from_df2_names.count(name) for name in set(missing_from_df2_names)}
    
    # Preparing data for plotting
    residues = list(set(exact_matches_names + missing_from_df1_names + missing_from_df2_names))
    hits_data = [hits_counts.get(res, 0) for res in residues]
    missing_df1_data = [-missing_df1_counts.get(res, 0) for res in residues]
    
    fig, ax = plt.subplots()

    # Plotting hits
    ax.bar(residues, hits_data, color='green', label='Hits')
    
    # Plotting misses from DF1
    ax.bar(residues, missing_df1_data, color='red', label='Missing from DF1')
    
    ax.set_ylabel('Counts')
    ax.set_title('Hits and Misses per Residue')
    ax.legend()

    plt.xticks(rotation=45)
    plt.tight_layout()  # Adjust layout to not cut off labels
    plt.show()

def print_summary(summary):
    """
    Print the summary of residue occurrences in a structured format for easy validation.
    Parameters:
        summary (dict): The summary dictionary returned by countingResiduesOccurrences function.
    """
    for category, residues_dict in summary.items():
        print(f"Category: {category}")
        for residue_name, details in residues_dict.items():
            print(f"  Residue Name: {residue_name}")
            for detail in details:
                if category == 'exact_matches' or category == 'missing_residues':
                    chain_id, residue_id = detail
                    print(f"    Chain ID: {chain_id}, Residue ID: {residue_id}")
                elif category == 'close_but_incorrect_matches' or category == 'incorrectly_present_residues':
                    chain_id, residue_id, other_residue = detail
                    print(f"    Chain ID: {chain_id}, Residue ID: {residue_id}, Other Residue: {other_residue}")
        print("\n")



def summary_to_dataframe(summary):
    """
    Convert the summary of residue occurrences into a pandas DataFrame with additional columns for missing residues.
    Parameters:
        summary (dict): The summary dictionary from countingResiduesOccurrences function.
    Returns:
        pd.DataFrame: A DataFrame with detailed residue information.
    """
    rows_list = []

    # Initialize counts for all residue names from exact matches, close but incorrect, and missing in df1
    all_residues = set(summary['exact_matches'].keys()) | set(summary['close_but_incorrect_matches'].keys()) | set(summary['missing_residues'].keys())

    for residue_name in all_residues:
        row_data = {
            'Residue Name': residue_name, 
            'ChainID': None,  # Placeholder, will update below
            'Hits': 0, 
            'Close But Incorrect Matches': 0, 
            'Incorrectly Present Residues': 0,
            'Missing Residues (from df1)': 0,
            'Missing Residues (from df2)': 0
        }

        # Update hits
        if residue_name in summary['exact_matches']:
            details = summary['exact_matches'][residue_name]
            row_data['Hits'] = len(details)
            row_data['ChainID'] = details[0][0]  # Assuming first detail's ChainID is representative

        # Update close but incorrect matches
        if residue_name in summary['close_but_incorrect_matches']:
            details = summary['close_but_incorrect_matches'][residue_name]
            row_data['Close But Incorrect Matches'] = len(details)
            if not row_data['ChainID']:
                row_data['ChainID'] = details[0][0]  # Update ChainID if not set

        # Update missing residues from df1
        if residue_name in summary['missing_residues']:
            details = summary['missing_residues'][residue_name]
            row_data['Missing Residues (from df1)'] = len(details)
            if not row_data['ChainID']:
                row_data['ChainID'] = details[0][0]  # Update ChainID if not set

        rows_list.append(row_data)

    # Handle missing residues from df2 separately
    for residue_name, details in summary['incorrectly_present_residues'].items():
        # Find or create the row for this residue
        for row in rows_list:
            if row['Residue Name'] == residue_name:
                row['Incorrectly Present Residues'] = len(details)
                row['Missing Residues (from df2)'] = len(details)  # Assuming same count as incorrectly present
                break
        else:
            # If not found, create a new row for this residue
            new_row = {
                'Residue Name': residue_name, 
                'ChainID': details[0][1],  # Using ChainID from incorrectly present details
                'Hits': 0, 
                'Close But Incorrect Matches': 0, 
                'Incorrectly Present Residues': len(details),
                'Missing Residues (from df1)': 0,
                'Missing Residues (from df2)': len(details)
            }
            rows_list.append(new_row)

    # Convert the list of rows into a DataFrame
    df_summary = pd.DataFrame(rows_list)

    # Fill missing ChainIDs with a placeholder if necessary (e.g., 'Unknown')
    df_summary['ChainID'].fillna('Unknown', inplace=True)

    return df_summary



def countingResiduesOccurrences2(df, residuesListGT, thresholdDistance, test_file):
    """
    Compute and return a DataFrame with precision and recall statistics for each type of residue.
    The computation is based on a DataFrame with residue correspondences and a threshold distance.
    """
    print('residues occurrences:')
    #thresholdDistance=5
    df_clean = df[df['distance'] <= thresholdDistance]
    print ('occupancy=', len(df_clean), ' of ', len(df), ' =>  ', format(len(df_clean)/len(df), '.2f') ,'%')
    print( 'compute for each residue:' )
    unique_residues = sorted(df_clean['residue1'].unique())
    counts = df_clean.groupby('residue1').apply(count_matches,'residue1','residue2')
    #print (counts)
    #print ('######')
    TP= (df_clean['residue1'] == df_clean['residue2']).sum()
    FP= (df_clean['residue1'] != df_clean['residue2']).sum()
    print ('TP=',TP, '  FP=', FP)

    df_stats = pd.DataFrame(columns=['residue', 'precision', 'recall'])
    for residue in unique_residues:
        #print ('#####  \n', residue)
        hits_per_residue_TP=counts[residue]['match']
        #founds_per_residue=sum(1 for item in residuesList1 if item[0] == residue)
        wrongly_assigned_FP= ((df['residue2'] == residue) & (~df['residue1'].isin([residue]))).sum()
        expected_per_residue=sum(1 for item in residuesListGT if item[0] == residue)
        miss_per_residue_FN=expected_per_residue-hits_per_residue_TP

        precision=hits_per_residue_TP/(hits_per_residue_TP+wrongly_assigned_FP) #TP/(TP+FP)
        recall=hits_per_residue_TP/(hits_per_residue_TP+miss_per_residue_FN) #TP/(TP+FN)

        print ('   ', residue , \
            '    HIT (TP) = ', "{:6}".format( hits_per_residue_TP ), \
            '    wrongly Assigned (FP) = ', "{:6}".format( wrongly_assigned_FP ) , \
            '    miss per residue (FN) =', "{:6}".format( miss_per_residue_FN ), \
            '    precision=', "{:6.2}".format( precision), \
            '    recall=', "{:6.2}".format( recall )    )
        result_row = {'residue': residue, 'precision': precision, 'recall': recall, 'test_file': test_file}
        df_stats = pd.concat([df_stats, pd.DataFrame(result_row, index=[0])], ignore_index=True)

    print( '######' )
    return df_stats


# function to plot stats
def plot_MultipleResiduesStats(df_stats, listNames, y_min=None, y_max=None):
    """
    Plot boxplots for the precision and recall values of multiple residue statistics.
    The input is a list of DataFrames, each containing precision and recall values for a type of residue, and a list of corresponding names.
    If y_min and y_max are provided, these are used as the limits for the y-axis. If not, the limits are set based on the data.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    cdf_precision = pd.DataFrame()
    cdf_recall = pd.DataFrame()
    sns.set_style("darkgrid")
    for idx in range (0, len(df_stats)):
        print (listNames[idx])
        label=listNames[idx]
        cdf_precision[label]=df_stats[idx]['precision']
        cdf_recall[label]=df_stats[idx]['recall']

    sns.boxplot(data=cdf_precision, ax=axs[0])
    sns.boxplot(data=cdf_recall, ax=axs[1])
    axs[0].set(title='Precision', ylabel='')
    axs[1].set(title='Recall', ylabel='')

    if y_min is None or y_max is None:
        y_min_0, y_max_0 = axs[0].get_ylim()
        y_min_1, y_max_1 = axs[1].get_ylim()
        y_min = min(y_min_0, y_min_1) if y_min is None else y_min
        y_max = max(y_max_0, y_max_1) if y_max is None else y_max

    axs[0].set_ylim(y_min, y_max)
    axs[1].set_ylim(y_min, y_max)

    fig.subplots_adjust(left=0.1)  # Adjust left margin to 14% of figure width

    plt.show()







def countingBackboneOccurrencies(df, atomListGT, thresholdDistance):
    """
    Compute and return a DataFrame with precision and recall statistics for each type of backbone atom.
    The computation is based on a DataFrame with backbone atom correspondences and a threshold distance.
    """
    print('atom occurrences:')
    #thresholdDistance=5
    # Filter the DataFrame based on distance
    df_clean = df[df['distance'] < thresholdDistance]
    #print (df_clean)

#        if atomList1[ii][0]=="CA" or atomList1[ii][0]=="N" or atomList1[ii][0]=="C":

    df_stats = pd.DataFrame(columns=['atom', 'precision', 'recall'])
    unique_atoms=['CA','C','N']
    for atom in unique_atoms:
        expected_atoms=sum(1 for item in atomListGT if item[0] == atom)
        counts = df_clean.groupby('atom1').apply(count_matches,'atom1','atom2')
        hits_atoms_TP=counts[atom]['match']
        wrongly_assigned_FP= ((df['atom2'] == atom) & (~df['atom1'].isin([atom]))).sum()
        miss_atom_FN=expected_atoms-hits_atoms_TP

        precision=hits_atoms_TP/(hits_atoms_TP+wrongly_assigned_FP) #TP/(TP+FP)
        recall=hits_atoms_TP/(hits_atoms_TP+miss_atom_FN) #TP/(TP+FN)

        print ('  atom: ', atom , \
            '    HIT (TP) = ', "{:6}".format( hits_atoms_TP ), \
            '    wrongly Assigned (FP) = ', "{:6}".format( wrongly_assigned_FP ) , \
            '    miss per residue (FN) =', "{:6}".format( miss_atom_FN ), \
            '    precision=', "{:6.2}".format( precision), \
            '    recall=', "{:6.2}".format( recall )    )
        result_row = {'atom': atom, 'precision': precision, 'recall': recall}
        df_stats = pd.concat([df_stats, pd.DataFrame(result_row, index=[0])], ignore_index=True)
    print( '######' )
    return df_stats



def plot_statistics(results_residues_list, inputTestFiles):
    # Concatenate all DataFrames into one
    df_stats_all = pd.concat(results_residues_list, ignore_index=True)

    # Define hydrophilicity based on the Kyte-Doolittle scale
    hydrophilicity = {'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
                      'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
                      'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2}

    # Define molecular weights and hydrophilicity
    molecular_weights = {'ALA': 89, 'ARG': 174, 'ASN': 132, 'ASP': 133, 'CYS': 121, 'GLN': 146, 'GLU': 147, 'GLY': 75,
                         'HIS': 155, 'ILE': 131, 'LEU': 131, 'LYS': 146, 'MET': 149, 'PHE': 165, 'PRO': 115, 'SER': 105,
                         'THR': 119, 'TRP': 204, 'TYR': 181, 'VAL': 117}

    hydrophilicity = {'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
                      'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
                      'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2}

    chou_fasman_parameters = {'ALA': 1.42, 'ARG': 0.94, 'ASN': 1.57, 'ASP': 1.60, 'CYS': 0.70, 'GLN': 1.56, 'GLU': 1.53,
                              'GLY': 0.57, 'HIS': 0.87, 'ILE': 1.08, 'LEU': 1.21, 'LYS': 0.74, 'MET': 1.45, 'PHE': 1.13,
                              'PRO': 0.55, 'SER': 1.41, 'THR': 1.25, 'TRP': 0.79, 'TYR': 1.25, 'VAL': 1.70}



    # Use this to sort residues
    residue_order = sorted(hydrophilicity.keys(), key=lambda x: hydrophilicity[x])

    # Plot
    fig, ax = plt.subplots(2,1, figsize=(10,12), sharex=True) # smaller size
    bar_width = 0.15  # Adjust as needed
    space_between_groups = 0.05  # Adjust as needed
    num_test_files = len(inputTestFiles)
    sns.set_palette("husl", num_test_files)

    for i, test_file in enumerate(inputTestFiles):
        precision_data = df_stats_all[(df_stats_all['test_file'] == test_file)]
        
        # Ensure all residues are present in the data
        for residue in residue_order:
            if residue not in precision_data['residue'].values:
                precision_data = precision_data.append({'residue': residue, 'precision': 0, 'recall': 0, 'test_file': test_file}, ignore_index=True)

        precision_data = precision_data.set_index('residue').loc[residue_order].reset_index()  # Sort the data

        # Adjust the x-values
        x_values = np.arange(len(precision_data['residue'])) + (i * (bar_width + space_between_groups))

        # Bar plot for precision
        ax[0].bar(x_values, precision_data['precision'], width=bar_width, label=f'Precision {test_file}')

        # Bar plot for recall
        ax[1].bar(x_values, precision_data['recall'], width=bar_width, label=f'Recall {test_file}')

    # Set common x-ticks and title
    ax[0].set_title('Precision per Residue for each Test File')
    ax[0].set_xticks(np.arange(len(precision_data['residue'])) + num_test_files * (bar_width + space_between_groups) / 2)
    ax[0].set_xticklabels(precision_data['residue'])
    ax[0].legend()

    ax[1].set_title('Recall per Residue for each Test File')
    ax[1].set_xticks(np.arange(len(precision_data['residue'])) + num_test_files * (bar_width + space_between_groups) / 2)
    ax[1].set_xticklabels(precision_data['residue'])
    ax[1].legend()

    plt.tight_layout()
    plt.show()





def main():

    #angelo_selection_best150k_14.pdb         angelo_selection_best150k_refined_25.pdb angelo_selection_full_9.pdb              pdb6g79.pdb
    parser = argparse.ArgumentParser(description='Compare two aligned pdb files to check how residues and backbone atoms are corresponding each others')
    parser.add_argument('--i', nargs='+', help='input test files')
    parser.add_argument('--labels', nargs='+', help='label name for each input file number')
    parser.add_argument('--gt', type=str, help='ground truth reference file')
    parser.add_argument('--mask', type=str, help='grayscale mask')
    parser.add_argument('--plot', action='store_true', help='Flag to enable plotting')
    parser.add_argument('--backbone', action='store_true', help='Analysis for Ca,C,N atoms')
    parser.add_argument('--radius', type=float, default=5.0, help='radius of search for a specific atom (default: 5.0A)')
    parser.add_argument('--ymin', type=float, default=None, help='Minimum y-value for the plot')  # new argument for minimum y-axis limit
    parser.add_argument('--ymax', type=float, default=None, help='Maximum y-value for the plot')  # new argument for maximum y-axis limit


    parser.set_defaults(plot=False)
    parser.set_defaults(backbone=False)

    args = parser.parse_args()
    
    if len(sys.argv) < 2:
        parser.print_help()
        exit(1)
    
    labels=[]

    inputTestFiles = args.i
    print( '========================' )
    print('Comparing two aligned pdb files to check how residues and/or backbone atoms are corresponding each others')
    for ii in range (0, len(inputTestFiles)):
        if labels and args.labels and len(args.labels) > ii:
            labels.append(args.labels[ii])
        else:
            labels.append(os.path.basename(inputTestFiles[ii][:inputTestFiles[ii].rfind('.')]))
        print('test file ',ii+1,'=      ', inputTestFiles[ii], '  label=', labels[ii])

    print('reference file = ', args.gt)
    print( '------------------------' )

    print('reference mask = ', args.mask)
    print( '------------------------' )
    #exit()

    residuesListGT=getResiduesList(args.gt)
    #exit()
    if (args.backbone):
        atomListGT=getBackboneAtomsList(args.gt)

    results_residues_list=[]
    results_atoms_list=[]
    for ii in range (0, len(inputTestFiles)):
        print( '========================' )
        print('processing test file ',ii+1,'=      ', inputTestFiles[ii])
        residuesList1=getResiduesList(inputTestFiles[ii])
        df1=getResiduesCorrespondenceMatrix(residuesListGT, residuesList1)
        df2=getResiduesCorrespondenceMatrix(residuesList1, residuesListGT)
        results=identify_mismatches_and_missing(df2, df1, thresholdDistance=3.5)
        plot_hits_and_misses(results)
        print (summary)
        exit
        #df_summary=summary_to_dataframe(summary)
        #print (df_summary)
        #print_summary(summary)
        #df_stats=countingResiduesOccurrences(df, residuesListGT, args.radius, inputTestFiles[ii])

        exit()
        df_stats['test_file'] = inputTestFiles[ii]
        results_residues_list.append(df_stats)
        if (args.backbone):
            atomList1=getBackboneAtomsList(inputTestFiles[ii])
            df_backbone=getBackboneCorrespondenceMatrix(atomList1, atomListGT)
            results_atoms_list=countingBackboneOccurrencies(df_backbone, atomListGT, args.radius)
        print( '------------------------' )

    if args.mask:
        maskPdbSelect(args.gt, 'tmp.pdb',  args.mask)

    
#    exit()
    if (args.plot):
        plot_MultipleResiduesStats(results_residues_list, labels, args.ymin, args.ymax)
        plot_statistics(results_residues_list, inputTestFiles)
    if args.plot and args.backbone:
        plot_MultipleResiduesStats(results_atoms_list, labels, args.ymin, args.ymax)




if __name__ == '__main__':
    main()
    

