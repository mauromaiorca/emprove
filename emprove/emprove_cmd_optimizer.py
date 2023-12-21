#!/usr/bin/python3


import argparse
import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns

emprove_parser = argparse.ArgumentParser(
    prog="emprove_optimizer",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
command = emprove_parser.add_subparsers(dest="command")





#################################
#################################
## predict_min_particles
def predict_min_particles(file_path="", outputImageFile="", showPlot=True, predicted_particles=None, ax=None):
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

        ax.set_title("Optimal Particle Number Estimation")
        ax.set_xlabel("Number of Particles")
        ax.set_ylabel("Mean Local Resolution Estimation Value")

        # Save to file if outputImageFile is provided
        if outputImageFile:
            directory = os.path.dirname(outputImageFile)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            plt.savefig(outputImageFile, dpi=300, format='png', bbox_inches='tight')
        
        # Display the plot
        if showPlot and ax is None:
            plt.show()

        if ax is None:
            plt.close()  # Close the plot to free memory

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
emprove_getNumParticles.add_argument("--plotOnFile", required=False, default="", type=str, help="Save the plot on file")
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
    result=predict_min_particles(args.locres, outputImageFile=args.plotOnFile, showPlot=args.plot)
    print(result)
    if not args.save=="":
        directory = os.path.dirname(args.save)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(args.save, 'w') as f:
            f.write(str(int(result)))


#################################
#################################
## automaticParticleSubsets Routine
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def skewed_gaussian(mean, skew_factor, size, seed=0, sampling_density_factor=0.15):
    #note: large scale translate into denser sampling
    np.random.seed(seed)  # Set the random seed
    values = np.random.normal(mean, mean*sampling_density_factor, size)
    values = mean + skew_factor * (values - mean)
    return values

def generate_particle_samples(predicted_particle_number, n_samples, min_val, max_val,sampling_density_factor=0.15):
    tmp_n_sample=n_samples
    final_samples=[]
    counter=0
    while len(final_samples) <= n_samples :
        counter+=1
        particle_samples = skewed_gaussian(predicted_particle_number, -0.7, tmp_n_sample*2, sampling_density_factor=sampling_density_factor)  # Generate extra to choose from
        cap_min=max(0, min_val)
        particle_samples = np.clip(particle_samples, cap_min, max_val)
        particle_samples = np.append(particle_samples, max_val)

        # Sort the particles
        particle_samples = sorted(particle_samples)

        # Ensure minimum separation
        min_gap = (predicted_particle_number - min_val) / (tmp_n_sample-2)
        final_samples = [int(particle_samples[-1])]  # Start with max value
        last_added = particle_samples[-1]
    
        for p in reversed(particle_samples):
            if abs(p - last_added) > min_gap and p not in final_samples:
            	final_samples.append(int(p))
            	last_added = p
            if len(final_samples) == tmp_n_sample:
            	break
        tmp_n_sample+=1

    return sorted(final_samples)


def automaticParticleSubsets2(locresResultsCsvFile, maxNumberOfParticles, number_of_sampling):
    if not locresResultsCsvFile == "":
    	df = pd.read_csv(locresResultsCsvFile)
    	predicted_particle_number = predict_min_particles(locresResultsCsvFile,showPlot=False)
    	min_val = df['numParticles'].min() * 2.0/3.0
    	max_val = df['numParticles'].max()
    	sampling_density_factor=0.15
    else:
    	predicted_particle_number = maxNumberOfParticles
    	min_val = maxNumberOfParticles * 1.0/3.0
    	max_val = maxNumberOfParticles
    	sampling_density_factor=0.25

    predicted_particles = generate_particle_samples(predicted_particle_number, number_of_sampling, min_val, max_val,sampling_density_factor=sampling_density_factor)
    return predicted_particles





#################################
## automaticParticleSubsets 
emprove_automaticParticleSubsets = command.add_parser (
    "automaticParticleSubsets", description="compute automaticParticleSubsets", help='compute automaticParticleSubsets'
)
emprove_automaticParticleSubsets.add_argument("--starFile", required=True, type=str, help="file with the input star file")
emprove_automaticParticleSubsets.add_argument("--locres", required=False, type=str, default="", help="file with the previous locres evaluation file is")
emprove_automaticParticleSubsets.add_argument("--save", required=False, type=str, default="", help="file where to save the particle to check in a csv fashon")
emprove_automaticParticleSubsets.add_argument("--plot", action="store_false", help="Display the plot")
emprove_automaticParticleSubsets.add_argument("--plotOnFile", required=False, default="", type=str, help="Save the plot on an image file")
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

    result=automaticParticleSubsets2(args.locres, num_non_null_items, 10)
    print ("result automaticParticleSubsets2=",result)
    exit(0)
    predict_min_particles(args.locres, showPlot=False,predicted_particles=result, outputImageFile=args.plotOnFile)
    print("Particles for performing selection",result)
    if not args.save == "":
    	directory = os.path.dirname(args.save)
    	if directory and not os.path.exists(directory):
    		os.makedirs(directory)
    	with open(args.save, 'w') as f:
    		f.write(','.join(map(str, result)))

    
def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "getNumParticles":
        getNumParticles(args)
    elif args.command == "automaticParticleSubsets" :
        automaticParticleSubsets(args)
    elif args.command == "logAnalyzer" :
        logAnalyzer(args)
    else:
        emprove_parser.print_help()


if __name__ == "__main__":
    main()

