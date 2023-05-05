import constants
from tools import models, plot

if __name__ == '__main__':
    params = constants.LIVING_DEATHRATE_PARAMS
    model = models.LinearModel(params)
    det_model = models.LinearModel(params)
    
    
    
    #l, v = det_model.get_stable_state()
    simulations_per_type = 20
    
    true_extinction_rates = model.calculate_extinction_rates().x
    simulated_extinction_rates = model.sample_extinction(simulations_per_type)
    #print(f'average growth rate {l} with stable proportions {v}')
    #print(f'Monte Carlo simulated extinction rates: {simulated_extinction_rates}')
    #print(f'How close (holistic) to a zero: {distance_from_zero}')
    #print(f'True extinction probability: {true_extinction_rates}')
    plot.plot_extinction_monte_carlo_graph(simulated_extinction_rates, true_extinction_rates, 
                                   simulations_per_type=simulations_per_type)