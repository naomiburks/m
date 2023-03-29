import constants
import plot
import models

if __name__ == '__main__':
    params = constants.TOY_PARAMS
    M = 5
    params[0] = M
    n_initial = constants.TOY_INITIAL(M)
    model = models.NoncollaborativeStochastic(params)
    data, timepoints = model.generate_timepoint_data(n_initial, 1, 100)
    plot.generate_graph(data, timepoints)