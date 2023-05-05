import constants
from tools import models, plot



if __name__ == '__main__':
   
    params = constants.LIVING_BIRTHRATE_PARAMS
    
    M = params['M']
    n_initial = constants.ALL_EQUAL(M, cell_count=500)
    model = models.LinearModel(params)
    
    data, timepoints = model.generate_timepoint_data(n_initial, 5, 200, stochastic=True, normalized=False)
    quasistable, growth = model.calculate_quasistable_distribution()

    plot.plot_simulation(data, timepoints, str(params))