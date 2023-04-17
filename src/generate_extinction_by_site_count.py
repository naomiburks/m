import constants
import models
import plot 


if __name__ == '__main__':
    base_params = constants.TOY_BARELY_LIVING_PARAMS
    data = []
    for M in [3, 10, 30, 100]:
        params = [param for param in base_params]
        params[0] = M
        model = models.NoncollaborativeStochastic(params)
        det_model = models.NoncollaborativeDeterministic(params)
        l, v = det_model.get_stable_state()
        true_extinction_rates = model.find_extinction_probabilities().x
        data.append(true_extinction_rates)
    other_params = base_params[1:]
    param_info = ('Other parameters:\n' + 
                f'b: {base_params[1]}\nd_0: {base_params[2]}\n' + 
                f'd_M: {base_params[3]}\nr_mu: {base_params[4]}\n' + 
                f'r_um: {base_params[5]}\np: {base_params[6]}')
    
    model = models.NoncollaborativeStochastic(base_params)
    limiting_probabilities = model.find_limit_extinction_probabilities(point_count = 1000)
    data.append(limiting_probabilities)
    plot.plot_extinction_limit_graph(data, param_info)