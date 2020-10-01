from model_driver import driver
from iotools import read_results
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # runs model
    outputfile = driver(create_ncf=True)
    print("Results",outputfile)
    # reads results from .nc-file
    results = read_results(outputfile)
    print("Type results",type(results))
    # plots ground water level for first ten nodes
    plt.figure()
    results['soil_ground_water_level'][:,0,:10].plot.line(x='date')
