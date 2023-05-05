import constants
from tools import models, plot


if __name__ == '__main__':
    #base_params = constants.LIVING_DEATHRATE_PARAMS
    base_params = constants.COLLABORATIVE_PARAMS
    data = []
    for M in [3, 10, 30, 100]:
        params = base_params.copy()
        params['M'] = M
        #model = models.ThresholdModel(params)
        model = models.SuperCollaborativeModel(params)
        det_model = models.NoncollaborativeDeterministic(params)
        l, v = det_model.get_stable_state()
        true_extinction_rates = model.calculate_extinction_rates().x
        data.append(true_extinction_rates)
    param_info = ('Other parameters:\n' + 
                f"b_0: {base_params['b_0']}\nb_M: {base_params['b_M']}\nd_0: {base_params['d_0']}\n" + 
                f"d_M: {base_params['d_M']}\nr_mu: {base_params['r_mu']}\n" + 
                f"r_um: {base_params['r_um']}\np: {base_params['p']}")
    
    model = models.SuperCollaborativeModel(base_params)
    limiting_probabilities = model.calculate_limit_extinction_rates(point_count = 1000)
    data.append(limiting_probabilities)
    plot.plot_extinction_limit_graph(data, param_info)