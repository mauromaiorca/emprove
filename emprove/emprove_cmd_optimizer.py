#!/usr/bin/python3


import argparse
import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns
import toml
import csv
import time
import shutil

emprove_parser = argparse.ArgumentParser(
    prog="emprove_optimizer",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
command = emprove_parser.add_subparsers(dest="command")





#################################
#################################
## predict_min_particles
def predict_min_particles(file_path="", outputImageFile="", outputSplineFile="", showPlot=True, predicted_particles=None, ax=None):
    sns.set_style("whitegrid")  # Set Seaborn's whitegrid style
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # Define brighter aesthetics using Seaborn's 'colorblind' palette
    colors = sns.color_palette("colorblind", 10)
    data_color = colors[0]    # Blue
    outlier_color = colors[2] # Orange
    spline_color = colors[3]  # Teal

    if file_path:
        df = pd.read_csv(file_path)
        x = df['numParticles'].values
        y = df['mean'].values

        # Detect outliers using IQR
        Q1 = np.percentile(y, 10)
        Q3 = np.percentile(y, 70)
        IQR = Q3 - Q1
        outlier_threshold_upper = Q3 + 1.5 * IQR
        outlier_threshold_lower = Q1 - 1.5 * IQR

        # Separate x and y values into fit and outlier sets
        fit_mask = (y >= outlier_threshold_lower) & (y <= outlier_threshold_upper)
        outlier_mask = ~fit_mask

        x_fit = x[fit_mask]
        y_fit = y[fit_mask]
        x_outlier = x[outlier_mask]
        y_outlier = y[outlier_mask]

        # Check if we have enough unique data points
        if len(x_fit) <= 1:
            print("Insufficient unique data points for spline fitting.")
            return None  # Or handle this scenario appropriately

        # Fit a spline with weights given to lower y-values
        k = min(4, len(x_fit) - 1)
        k = max(1, min(4, len(x_fit) - 1))
        weights = 1 / y_fit
        s = UnivariateSpline(x_fit, y_fit, w=weights, k=k, s=len(x_fit))
        x_smooth = np.linspace(min(x_fit), max(x_fit), 1000)
        y_smooth = s(x_smooth)

        # Create whisker plots for fitting data points
        whisker_data = df[['max', 'highQuartile', 'mean', 'lowQuartile', 'min']].values.T
        ax.errorbar(x_fit, y_fit, yerr=[y_fit - whisker_data[4, fit_mask], whisker_data[0, fit_mask] - y_fit],
                    fmt='o', color=data_color, label='Local Resolution Variation', capsize=5, capthick=1.5)

        # Create whisker plots for outliers
        ax.errorbar(x_outlier, y_outlier, yerr=[y_outlier - whisker_data[4, outlier_mask], whisker_data[0, outlier_mask] - y_outlier],
                    fmt='o', color=outlier_color, label='Outliers', capsize=5, capthick=1.5)

        # Plot the fitted spline
        ax.plot(x_smooth, y_smooth, label='Fitted Spline', color=spline_color, linewidth=2.5)

        # Determine the x-value corresponding to the minimum y-value on the spline
        idx = np.argmin(y_smooth)
        predicted_particle_number = x_smooth[idx]

        # Plot a vertical line where the predicted_particle_number is with a label for the legend
        ax.axvline(x=predicted_particle_number, color=colors[1], linestyle='--', linewidth=1.5, label="Estimation")  # A different shade of colorblind palette for estimation line

        # Add the estimated number as a top x-tick
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks([predicted_particle_number])
        ax2.set_xticklabels([f'{int(predicted_particle_number)}'])

        # If the predicted_particle_list is provided, plot vertical dotted lines for each prediction
        if predicted_particles:
            for particle in predicted_particles:
                ax.axvline(x=particle, color=colors[4], linestyle=':', linewidth=0.8)  # Another color from the palette

        # Adjust x-axis limits
        if predicted_particles:
            ax.set_xlim(min(min(x), min(predicted_particles)), max(max(x), max(predicted_particles)))

        # Floating legend with opaque background
        leg = ax.legend(loc='best', frameon=True)
        leg.get_frame().set_alpha(1.0)  # Non-transparent legend

        ax.set_title("Mean Local Resolution per Number of Particles")
        ax.set_xlabel("Number of Particles")
        ax.set_ylabel("Mean Local Resolution Estimation Value")

        # Save to file if outputImageFile is provided
        if outputImageFile:
            directory = os.path.dirname(outputImageFile)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            #print("saving on file ",outputImageFile)
            plt.savefig(outputImageFile, dpi=300, format='png', bbox_inches='tight')
        
        # Display the plot
        if showPlot:
            plt.show()
        plt.close()  # Close the plot to free memory

        if outputSplineFile:
            spline_data = pd.DataFrame({
                'numParticles': x_smooth,
                'estimatedMeanResolution': y_smooth
            })
            # Ensure the output directory exists
            output_directory = os.path.dirname(outputSplineFile)
            if output_directory and not os.path.exists(output_directory):
                os.makedirs(output_directory)
            spline_data.to_csv(outputSplineFile, index=False)   

        return int(predicted_particle_number) if file_path else None
    return None



def create_summary_figure(basename, rows=3, columns=3):
    import matplotlib.gridspec as gridspec
    print("basename=",os.path.basename(basename))
    def extract_number(dir_name, prefix=os.path.basename(basename)):
        # Get the part of the string after the known prefix
        if prefix in dir_name:
            part_after_prefix = dir_name.split(prefix)[-1]  # This gets the part of the string after the prefix
            # Now, extract consecutive digits from the beginning of this substring
            number_str = ""
            for char in part_after_prefix:
                if char.isdigit():
                    number_str += char
                else:
                    break  # Break at the first non-digit character
            
            if number_str:
                return int(number_str)
        
        return 0  # Return 0 if no number is found or handle it differently

    base_folder = '.'
    if os.path.sep in basename:
        base_folder = os.path.dirname(basename)
    all_directories = [d for d in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, d)) and "_emprove_SCI__" in d]

    # Sort directories based on the embedded number
    all_directories.sort(key=extract_number)

    # Break the selected directories into chunks for multiple figures
    figure_count = 0
    for i in range(0, len(all_directories), rows * columns):
        figure_count += 1

        selected_dirs = all_directories[i:i+rows*columns]

        # Create a figure and GridSpec object
        fig = plt.figure(figsize=(16, 5 * rows))
        gs = gridspec.GridSpec(rows, columns)  # Grid based on input rows and columns

        for j, dir_name in enumerate(selected_dirs):
            file_path = os.path.join(base_folder, dir_name, "bestRanked_locres_values.csv")

            ax = fig.add_subplot(gs[j])

            # Call your plotting function using the given axes
            predict_min_particles(file_path=file_path, showPlot=False, ax=ax)
            iteration_number = extract_number(dir_name)
            ax.set_title(f"Iteration {iteration_number}")  # Use the extracted number as title for each subplot

        plt.tight_layout()  # Adjust layout

        # Save the figure with an incremented name for multiple figures
        fig.savefig(os.path.join(base_folder, f"summary_{figure_count}.png"), dpi=300, format='png', bbox_inches='tight')
        fig.savefig(os.path.join(base_folder, f"summary_{figure_count}.pdf"), format='pdf', bbox_inches='tight')
        plt.show()

#################################
## log analyzer
emprove_logAnalyzer = command.add_parser (
    "logAnalyzer", description="compute logAnalyzer", help='logAnalyzer'
)
emprove_logAnalyzer.add_argument("--prefix", required=True, type=str, help="folders prefix, e.g. _emprove_SCI__1.00_scored_selection_")
emprove_logAnalyzer.add_argument("--numRows", required=False, type=int, default=3, help="now of rows in reported figure (default=3)")
emprove_logAnalyzer.add_argument("--numColumns", required=False, type=int, default=3, help="now of columns in reported figure (default=3)")
def logAnalyzer(args):
    create_summary_figure(args.prefix, args.numRows, args.numColumns)

#################################
## Analyse reconstruction script
emprove_getNumParticles = command.add_parser (
    "getNumParticles", description="compute the optimal Target Number Of Particles, it requires the locres summary file", help='compute the optimal Target Number Of Particles, it requires the locres summary file'
)
emprove_getNumParticles.add_argument("--locres", required=True, type=str, help="file with locres evaluation")
emprove_getNumParticles.add_argument("--plot", action="store_true", help="Display the plot")
emprove_getNumParticles.add_argument("--mean_res", action="store_true", help="predict the best mean local resolution")
emprove_getNumParticles.add_argument("--plotOnFile", required=False, default="", type=str, help="Save the plot on file")
emprove_getNumParticles.add_argument("--saveSplineOnCsv", required=False, default="", type=str, help="Save the plot on file")
emprove_getNumParticles.add_argument("--save", required=False, default="", type=str, help="Save the best reconstruction particle on file")
def getNumParticles(args):
    
    if (not os.path.isfile(args.locres)):
        print('ERROR: file \"',args.locres,'\" not existing')
        exit()
    #df = pd.read_csv(args.locres)
    # Sorting the DataFrame
    #df_sorted = df.sort_values(by=['mean', 'numParticles', 'max', 'highQuartile', 'lowQuartile', 'min'], ascending=True)
    # Selecting the "numParticles" from the first row
    #result = int(df_sorted.iloc[0]['numParticles'])
    result=predict_min_particles(args.locres, outputImageFile=args.plotOnFile, outputSplineFile=args.plotOnFile, showPlot=args.plot)
    if (args.mean_res):
        data = pd.read_csv(args.locres)  # Replace 'your_file.csv' with your file path
        idx_min_mean = data['mean'].idxmin()
        result = data.loc[idx_min_mean, 'numParticles']
    print(result)

    if not args.save=="":
        directory = os.path.dirname(args.save)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(args.save, 'w') as f:
            f.write(str(int(result)))



#################################
## Analyse reconstruction script
emprove_plotOverview = command.add_parser (
    "plotOverview", description="plot the overview", help='plot the overview'
)
emprove_plotOverview.add_argument("--overview", required=True, type=str, help="overviewFile")
emprove_plotOverview.add_argument("--o", required=False, default="", type=str, help="output png file")
def plotOverview(args):
    if (not os.path.isfile(args.overview)):
        print('ERROR: overview file \"',args.overview,'\" not existing')
        exit()
    data=overview_read_item_at_iteration(args.overview, 1)
    print (data)



#################################
#################################
## automaticParticleSubsets Routine
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def skewed_gaussian(mean, skew_factor, size, seed=0, sampling_density_factor=0.15):
    #note: large scale translate into denser sampling
    #skew_factor=value that determines the direction and magnitude of skewness in the generated distribution. Positive values skew the distribution to the right, and negative values skew it to the left.
    np.random.seed(seed)
    values = np.random.normal(mean, mean * sampling_density_factor, size)
    values = mean + skew_factor * (values - mean)
    return values

def generate_particle_samples(predicted_particle_number, n_samples, min_val, max_val, sampling_density_factor=0.15, seed=0):
    final_samples = []
    counter = 0
    global_max_val = max_val
    skew_factor=-1.5
    max_val = (max_val + 2 * predicted_particle_number) / 3.0  # Adjust max_val for skewness

    while len(final_samples) < n_samples - 3:  # Adjust to leave space for the smallest and largest values
        counter += 1
        particle_samples = skewed_gaussian(predicted_particle_number, skew_factor, n_samples * 3, seed=seed, sampling_density_factor=sampling_density_factor)
        particle_samples = np.clip(particle_samples, max(0, min_val), max_val)

        # Sort and filter for minimum gap, avoiding duplicates, and preserving range
        particle_samples = np.sort(particle_samples)
        last_added = float('-inf')  # Initialize with a value that is always smaller
        for p in particle_samples:
            if len(final_samples) >= n_samples - 3 or p >= global_max_val:
                break
            if (p - last_added) >= (predicted_particle_number - min_val) / (n_samples - 1) and p not in final_samples:
                final_samples.append(p)
                last_added = p

    # Ensure to include the largest values and predicted_particle_number/2
    final_samples = [predicted_particle_number/2] + [predicted_particle_number]+final_samples + [global_max_val]

    # Convert to integers
    final_samples = [int(value) for value in final_samples]
    final_samples = sorted(set(final_samples))

    # Sort final samples and ensure we only have the desired number of samples, in case of duplicates near boundaries
    return final_samples



def automaticParticleSubsetsCore(locresResultsCsvFile, maxNumberOfParticles, number_of_sampling, randomSeed=True, showPlot=False, outputImageFile=""):
    if not locresResultsCsvFile == "":
        df = pd.read_csv(locresResultsCsvFile)
        predicted_particle_number = predict_min_particles(locresResultsCsvFile,  showPlot=showPlot, outputImageFile=outputImageFile)
        print("predicted_particle_number = ", predicted_particle_number )
        min_val = df['numParticles'].min() * 2.0/3.0
        max_val = df['numParticles'].max()
        sampling_density_factor=0.25
    else:
        predicted_particle_number = maxNumberOfParticles
        min_val = maxNumberOfParticles * 1.0/3.0
        max_val = maxNumberOfParticles
        sampling_density_factor=0.25
    if (randomSeed):
        seed = int(time.time())
    else:
        seed = 0
     
    predicted_particles = generate_particle_samples(predicted_particle_number, number_of_sampling, min_val, max_val,sampling_density_factor=sampling_density_factor, seed=seed)
    return predicted_particles





#################################
## automaticParticleSubsets 
emprove_automaticParticleSubsets = command.add_parser (
    "automaticParticleSubsets", description="compute automaticParticleSubsets", help='compute automaticParticleSubsets'
)
emprove_automaticParticleSubsets.add_argument("--starFile", required=True, type=str, help="file with the input star file")
emprove_automaticParticleSubsets.add_argument("--locres", required=False, type=str, default="", help="file with the previous locres evaluation file is")
emprove_automaticParticleSubsets.add_argument("--save", required=False, type=str, default="", help="file where to save the particle to reconstruct and compute locres, stored as comma separated file (csv)")
emprove_automaticParticleSubsets.add_argument("--plot", action="store_true", help="Display the plot")
emprove_automaticParticleSubsets.add_argument("--plotPrediction", action="store_true", help="Display the predictions for next plot")
emprove_automaticParticleSubsets.add_argument("--plotOnFile", required=False, default="", type=str, help="Save the plot on an image file")
emprove_automaticParticleSubsets.add_argument("--numSamples", required=False, type=int, default=10,  help="number of samples")
#emprove_automaticParticleSubsets.add_argument("--pureRandom", action="store_true", help="pure random number, not reproducible")
def automaticParticleSubsets(args):
    if (not os.path.isfile(args.starFile)):
        print('ERROR: file \"',args.starFile,'\" not existing')
        exit()
    from emprove import starHandler
    starFile=starHandler.readStar(args.starFile)

    num_non_null_items = int(starFile.count()[0])

    if os.path.isfile(args.locres):
        expectedEstimatedParticlesNumber=predict_min_particles(args.locres,showPlot=False)
        print("expected number of particles=",expectedEstimatedParticlesNumber)
    else:
        args.locres=""
        print ('WARNING: you might want to specify a valid locres file, however it is ok for the first iteration\n')
        expectedEstimatedParticlesNumber=num_non_null_items

    #print ("plot on file=",args.plotOnFile)
    result=automaticParticleSubsetsCore(args.locres, num_non_null_items, args.numSamples, showPlot=args.plot, outputImageFile=args.plotOnFile)
    print ("result automaticParticleSubsetsCore=",result)
#    predict_min_particles(args.locres, showPlot=args.plot, predicted_particles=result, outputImageFile=args.plotOnFile)
#    print("Particles for performing selection",result)
    if not args.save == "":
        directory = os.path.dirname(args.save)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(args.save, 'w') as f:
            f.write(','.join(map(str, result)))



#####################################################
######## OVERVIEW FILE
#####################################################
def overview_read_item_at_iteration(overview_file_path, iteration_number):
    # Load the TOML file
    try:
        data = toml.load(overview_file_path)
    except Exception as e:
        print(f"Error reading TOML file: {e}")
        return None

    # Iterate through each section to find the matching iteration
    for section in data:
        if "iteration" in data[section] and data[section]["iteration"] == iteration_number:
            return data[section].get("target_num_particles", None)

    print(f"Iteration {iteration_number} not found in the file.")
    return None



def remove_duplicates_toml(file_path):
    # Read the TOML file as plain text
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Process each line to remove duplicates, keeping only the last occurrence
    new_lines = []
    seen_keys = set()
    for line in reversed(lines):
        if '=' in line:
            key = line.split('=', 1)[0].strip()
            if key not in seen_keys:
                seen_keys.add(key)
                new_lines.append(line)
        else:
            new_lines.append(line)
            seen_keys.clear()  # Reset for a new table

    # Reverse the lines back to original order
    new_lines = new_lines[::-1]

    # Write the processed lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(new_lines)



emprove_generate_overview = command.add_parser (
    "generate_overview", description="generate overview", help='generate overview'
)
emprove_generate_overview.add_argument("--directory", required=True, type=str, help="directory with the overview to be generated")
emprove_generate_overview.add_argument("--verbose", action="store_true", help="verbose overview display")
def generate_overview(args):
    prefix = "_emprove_"
    # List all entries in the directory
    entries = os.listdir(args.directory)
    #print ("args.directory=",args.directory)

    # Filter directories that start with the prefix and end with a number
    emprove_dirs = []
    for entry in entries:
        if entry.startswith(prefix) and "_" in entry:
            parts = entry.split("_")
            if parts[-1].isdigit():  # Check if the last part is a number
                emprove_dirs.append(entry)

    # Sort the directories based on the numeric part
    emprove_dirs.sort(key=lambda x: int(x.split("_")[-1]))
    #print (emprove_dirs)

    fileSettings=os.path.join(args.directory,"session_settings.toml")
    with open(fileSettings, 'r') as file:
        dataSettings = toml.load(fileSettings)
    sigmaInitial=round(float(dataSettings.get("sigma", 1)),2)
    sigmaCurrent=sigmaInitial
    minimum_sigma_allowed=float(dataSettings.get("minimum_sigma_allowed", 0.6))
    sigma_decreasing_step=float(dataSettings.get("sigma_decreasing_step", 0.05))


    output_log_file_string=""
    lastSelectionsNoImprovements=0
    for index, ii in enumerate(emprove_dirs):
        if (args.verbose):
            print ("VERBOSE: checking directory ", os.path.join(args.directory, ii))
        bestRanked_path = os.path.join(args.directory, ii, "bestRanked_locres_values.csv")
        targetNum_path = os.path.join(args.directory, ii, "target_num_of_particles.csv")
        if not os.path.exists(bestRanked_path) or not os.path.exists(targetNum_path):
            print("WARNING: missing relevant files in the directory ",os.path.join(args.directory, ii),", ignoring it ")
            continue

        with open(bestRanked_path, mode='r') as file:
            reader = csv.DictReader(file)
            data1 = [row for row in reader]

        with open(targetNum_path, mode='r') as file:
            reader = csv.DictReader(file)
            data2 = [row for row in reader]

        # Combine, remove duplicates, and sort
        combined_data = {f"{row['numParticles']}_{row['mean']}": row for row in data1 + data2}.values()
        combined_data = sorted(combined_data, key=lambda x: int(x['numParticles']))

        #write combined data on a file
        fullPrediction_locres_file = os.path.join(args.directory,ii,"fullPrediction_bestRanked_locres_values.csv")
        headers = combined_data[0].keys() if combined_data else []
        with open(fullPrediction_locres_file, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=headers)
            writer.writeheader()
            writer.writerows(combined_data)



        # Find the row with the lowest mean
        lowest_mean = float('inf')
        lowest_mean_row = None
        for row in combined_data:
            if float(row['mean']) < lowest_mean:
                lowest_mean = float(row['mean'])
                lowest_mean_row = row

        # Find the row with the largest 'numParticles'
        largest_num_particles = 0
        largest_num_particles_row = None  # Define the variable
        for row in combined_data:
            if int(row['numParticles']) > largest_num_particles:
                largest_num_particles = int(row['numParticles'])
                largest_num_particles_row = row  # Assign the row data

        mean_of_lowest_mean = lowest_mean_row['mean']
        numParticles_of_lowest_mean = lowest_mean_row['numParticles']
        mean_of_largest_num_particles = largest_num_particles_row['mean']

        # Convert the dictionary rows to comma-separated strings
        lowest_mean_row_str = ','.join(map(str, lowest_mean_row.values()))
        largest_num_particles_str = ','.join(map(str, largest_num_particles_row.values()))

        # Print the results as comma-separated strings
        if (args.verbose):
            print("   lowest_mean_row:", lowest_mean_row_str)
            print("   mean_of_lowest_mean:", mean_of_lowest_mean)
            print("   largest_num_particles_row:", largest_num_particles_str)
            print("   mean_of_largest_num_particles:", mean_of_largest_num_particles)
            print("   numParticles_of_lowest_mean:", numParticles_of_lowest_mean)

        
        toCheckDir=os.path.join(args.directory,ii)
        targetNumOfParticles=0


        selection_block='tag="'+ii+'"\n'
        selection_block+='working_directory = "'+toCheckDir+'"\n'
        selection_block+="selection_number = "+str(int(ii.split("_")[-1]))+"\n"
        selection_block+="reference_num_particles = "+str(numParticles_of_lowest_mean)+"\n"
        selection_block+='reference_starFile = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'.star')+'"\n'
        selection_block+='reference_mapA = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'_recH1.mrc')+'"\n'
        selection_block+='reference_mapB = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'_recH2.mrc')+'"\n'
        selection_block+='reference_locres_stats = "'+lowest_mean_row_str+'"\n'
        selection_block+="reference_locres_mean = "+mean_of_lowest_mean+"\n"        
        selection_block+='computed_locres_file = "'+os.path.join(args.directory,ii,'fullPrediction_bestRanked_locres_values.csv')+'"\n'

        if index == 0:
            output_log_file_string+="\n[[_emprove_selection_0]]\n"
            output_log_file_string+='tag="'+ii+'"\n'
            output_log_file_string+='working_directory = "'+toCheckDir+'"\n'
            output_log_file_string+="reference_num_particles="+str(largest_num_particles)+"\n"
            output_log_file_string+='reference_starFile = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'.star')+'"\n'
            output_log_file_string+='reference_mapA = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'_recH1.mrc')+'"\n'
            output_log_file_string+='reference_mapB = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'_recH2.mrc')+'"\n'
            output_log_file_string+='reference_locres_stats = "'+largest_num_particles_str+'"\n'
            output_log_file_string+="reference_locres_mean = "+mean_of_largest_num_particles+"\n"
            output_log_file_string+='computed_locres_file = "'+os.path.join(args.directory,ii,'bestRanked_locres_values.csv')+'"\n'
            output_log_file_string+="improve_previouses_locres_mean = true\n"
            output_log_file_string+='session_settings = "'+os.path.join(args.directory,'session_settings.toml')+'"\n'
            #output_log_file_string+='SCI_sigma = "'+"{:.2f}".format(sigmaCurrent)+'"\n'
            #output_log_file_string+='SCI_sigma_next_selection = "'+"{:.2f}".format(sigmaCurrent)+'"\n'
            target_reference_starFile = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'.star')
            target_reference_mapA = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'_recH1.mrc')
            target_reference_mapB = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(largest_num_particles)+'_recH2.mrc')
            best_locres=mean_of_largest_num_particles
            best_selectionBlock=selection_block
            best_selectionBlock_size=largest_num_particles

        improve_previouses_locres_mean='false'
        if mean_of_lowest_mean < best_locres:
            best_locres=mean_of_lowest_mean
            best_selectionBlock=selection_block
            best_selectionBlock_size=numParticles_of_lowest_mean
            improve_previouses_locres_mean='true'
            target_reference_starFile = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'.star')
            target_reference_mapA = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'_recH1.mrc')
            target_reference_mapB = os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(numParticles_of_lowest_mean)+'_recH2.mrc')
            lastSelectionsNoImprovements=0
        else:
            improve_previouses_locres_mean='false'
            lastSelectionsNoImprovements+=1
            #if sigmaCurrent>minimum_sigma_allowed:
            #    sigmaCurrent=sigmaCurrent-sigma_decreasing_step
        #best_selectionBlock+='SCI_sigma = "'+"{:.2f}".format(sigmaCurrent)+'"\n'


        output_log_file_string+="\n[[_emprove_selection_"+str(int(ii.split("_")[-1]))+"]]\n"
        output_log_file_string+=selection_block
        #output_log_file_string+='SCI_sigma = "'+"{:.2f}".format(sigmaCurrent)+'"\n'
        output_log_file_string+="improve_previouses_locres_mean = "+improve_previouses_locres_mean+"\n"

        output_log_file_string+="\n"

    #shutil.copy(target_reference_starFile, os.path.join(args.directory))
    #shutil.copy(target_reference_mapA, os.path.join(args.directory))
    #shutil.copy(target_reference_mapB, os.path.join(args.directory))


    shutil.copy(target_reference_starFile, os.path.join(args.directory,'reference_subset.star'))
    shutil.copy(target_reference_mapA, os.path.join(args.directory,'reference_subset_mapA.mrc'))
    shutil.copy(target_reference_mapB, os.path.join(args.directory,'reference_subset_mapB.mrc'))

    percentage_particles_retained = 0
    if largest_num_particles > 0:
        percentage_particles_retained = (float(best_selectionBlock_size)/float(largest_num_particles))*100



    #Writing to a TOML file
    #    toml.dump(dir_data, toml_file)
    with open(os.path.join(args.directory,'overview.txt'), 'w') as toml_file:
        toml_file.write("#Selections overview\n\n")
        toml_file.write("[[_emprove_target_selection]]\n")
        toml_file.write(best_selectionBlock)
        toml_file.write("last_consecutive_non_improving_selections= "+str(lastSelectionsNoImprovements)+"\n")
        toml_file.write("percentage_particles_retained = "+"{:.2f}".format(float(percentage_particles_retained))+"\n")
        toml_file.write(output_log_file_string)
    
    remove_duplicates_toml(os.path.join(args.directory,'overview.txt'))
    return emprove_dirs







def generate_overview2(args):
    prefix = "_emprove_"
    # List all entries in the directory
    entries = os.listdir(args.directory)
    #print ("args.directory=",args.directory)

    # Filter directories that start with the prefix and end with a number
    emprove_dirs = []
    for entry in entries:
        if entry.startswith(prefix) and "_" in entry:
            parts = entry.split("_")
            if parts[-1].isdigit():  # Check if the last part is a number
                emprove_dirs.append(entry)

    # Sort the directories based on the numeric part
    emprove_dirs.sort(key=lambda x: int(x.split("_")[-1]))
    #print (emprove_dirs)

    output_log_file_string=""
    

    for index, ii in enumerate(emprove_dirs):
        if (args.verbose):
            print (ii)
        toCheckDir=os.path.join(args.directory,ii)
        targetNumOfParticles=0
        with open(os.path.join(args.directory,ii,"target_num_of_particles.csv"),'r') as file:
            first_line = file.readline()
            for word in first_line.split():
                if word.isdigit():
                    targetNumOfParticles=int(word)
        #read locres values for the target and put them in a file
        target_locres_values=" "
        target_locres_mean = "9999"
        global_num_particles=0
        global_locres_values=" "
        global_locres_mean = "9999"

        with open(os.path.join(args.directory,ii,'target_locres_values.csv'), mode='r') as tmp_locres_file:
            csv_reader = csv.reader(tmp_locres_file)
            headers = next(csv_reader)
            mean_index = headers.index('mean')
            target_all_values = []
            for values in csv_reader:
                target_all_values.append(values)
                target_locres_mean = values[mean_index]
            target_locres_values = ' '.join([','.join(row) for row in target_all_values])
        with open(os.path.join(args.directory,ii,'bestRanked_locres_values.csv'), mode='r') as tmp_locres_initial_file:
            csv_reader = csv.reader(tmp_locres_initial_file)
            headers = next(csv_reader)
            mean_index = headers.index('mean')
            num_particles_index = headers.index('numParticles')
            max_num_particles = 0
            max_row = None
            best_local_resolution = 9999
            best_local_resolution_row = None

            all_values = []
            for values in csv_reader:
                num_particles = int(values[num_particles_index])
                local_resolution = float(values[mean_index])
                if num_particles > max_num_particles:
                    max_num_particles = num_particles
                    max_row = values
                if local_resolution < best_local_resolution:
                    best_local_resolution = local_resolution
                    best_local_resolution_row = values                    
            if max_row:
                global_locres_mean = max_row[mean_index]
                global_num_particles = max_row[num_particles_index]
                global_locres_values = ','.join(max_row)


        selection_block='tag="'+ii+'"\n'
        selection_block+='working_directory = "'+toCheckDir+'"\n'
        selection_block+="selection_number = "+str(int(ii.split("_")[-1]))+"\n"
        selection_block+="reference_num_particles = "+str(targetNumOfParticles)+"\n"
        selection_block+='reference_starFile = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(targetNumOfParticles)+'.star')+'"\n'
        selection_block+='reference_mapA = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(targetNumOfParticles)+'_recH1.mrc')+'"\n'
        selection_block+='reference_mapB = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(targetNumOfParticles)+'_recH2.mrc')+'"\n'
        selection_block+='reference_locres_stats = "'+target_locres_values+'"\n'
        selection_block+="reference_locres_mean = "+target_locres_mean+"\n"        
        selection_block+='computed_locres_file = "'+os.path.join(args.directory,ii,'bestRanked_locres_values.csv')+'"\n'


        if index == 0:
            output_log_file_string+="\n[[_emprove_selection_0]]\n"
            output_log_file_string+='tag="'+ii+'"\n'
            output_log_file_string+='working_directory = "'+toCheckDir+'"\n'
            output_log_file_string+="reference_num_particles="+str(global_num_particles)+"\n"
            output_log_file_string+='reference_starFile = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(global_num_particles)+'.star')+'"\n'
            output_log_file_string+='reference_mapA = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(global_num_particles)+'_recH1.mrc')+'"\n'
            output_log_file_string+='reference_mapB = "'+os.path.join(args.directory,ii,'norm_'+ii+'_best'+str(global_num_particles)+'_recH2.mrc')+'"\n'
            output_log_file_string+='reference_locres_stats = "'+global_locres_values+'"\n'
            output_log_file_string+="reference_locres_mean = "+global_locres_mean+"\n"
            output_log_file_string+='computed_locres_file = "'+os.path.join(args.directory,ii,'bestRanked_locres_values.csv')+'"\n'
            output_log_file_string+="improve_previouses_locres_mean = true\n"
            best_locres=global_locres_mean
            best_selectionBlock=selection_block

        improve_previouses_locres_mean='false'
        if target_locres_mean < best_locres:
            best_locres=target_locres_mean
            best_selectionBlock=selection_block
            improve_previouses_locres_mean='true'



        output_log_file_string+="\n[[_emprove_selection_"+str(int(ii.split("_")[-1]))+"]]\n"
        output_log_file_string+=selection_block
        output_log_file_string+="improve_previouses_locres_mean = "+improve_previouses_locres_mean+"\n"
        output_log_file_string+="\n"


    #Writing to a TOML file
    #    toml.dump(dir_data, toml_file)
    with open(os.path.join(args.directory,'overview.txt'), 'w') as toml_file:
        toml_file.write("#Selections overview\n\n")
        toml_file.write("[[_emprove_target_selection]]\n")
        toml_file.write(best_selectionBlock)
        toml_file.write(output_log_file_string)
    return emprove_dirs




emprove_getTarget = command.add_parser (
    "getTarget", description="get the target from a specific overview file", help='get the target from a specific overview file'
)
emprove_getTarget.add_argument("--overviewFile", required=True, type=str, help="toml overview file to acquire target information")
emprove_getTarget.add_argument("--particles", action="store_true", help="gets the target star file with relevant particles")
emprove_getTarget.add_argument("--map1", action="store_true", help="gets the target map1 file")
emprove_getTarget.add_argument("--map2", action="store_true", help="gets the target map2 file")
emprove_getTarget.add_argument("--stats", action="store_true", help="num particles, min, quartile, mean, quartile, max")
#emprove_getTarget.add_argument("--sigma", action="store_true", help="get the sigma for the next iteration")
def getTarget(args):
    data = toml.load(args.overviewFile)
    target_selection = data["_emprove_target_selection"][0]  # Access the first element
    if args.particles:
        target_starfile = target_selection["reference_starFile"]
        print(target_starfile)
    if args.map1:
        reference_mapA = target_selection["reference_mapA"]
        print(reference_mapA)
    if args.map2:
        reference_mapB = target_selection["reference_mapB"]
        print(reference_mapB)
    if args.stats:
        reference_stats = target_selection["reference_locres_stats"]
        print(reference_stats)
    #if args.sigma:
    #    reference_sigma = target_selection["SCI_sigma"]
    #    print(reference_sigma)
    
def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "getNumParticles":
        getNumParticles(args)
    elif args.command == "automaticParticleSubsets" :
        automaticParticleSubsets(args)
    elif args.command == "logAnalyzer" :
        logAnalyzer(args)
    elif args.command == "plotOverview" :
        plotOverview(args)
    elif args.command == "generate_overview" :
        generate_overview(args)
    elif args.command == "getTarget" :
        getTarget(args)
    else:
        emprove_parser.print_help()


if __name__ == "__main__":
    main()

