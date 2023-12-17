import os.path

def getTargetNumberOfParticles(fileWithLocres):
    print("getTargetNumberOfParticles")
    exit()
    import pandas as pd

    if (not os.path.isfile(fileWithLocres)):
        print('ERROR: file \"',fileWithLocres,'\" not existing')
        exit()
    df = pd.read_csv(fileToAnalyse)
    # Sorting the DataFrame
    df_sorted = df.sort_values(by=['mean', 'numParticles', 'max', 'highQuartile', 'lowQuartile', 'min'], ascending=True)
    # Selecting the "numParticles" from the first row
    result = int(df_sorted.iloc[0]['numParticles'])
    print(result)
