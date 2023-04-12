import constants
import models
import plot 


if __name__ == '__main__':
    params = constants.TOY_LIVING_PARAMS
    model = models.NoncollaborativeStochastic(params)
    det_model = models.NoncollaborativeDeterministic(params)
    l, v = det_model.get_stable_state()
    simulations_per_type = 200
    simulated_extinction_rates = model.sample_extinction(simulations_per_type)
    true_extinction_rates = model.find_extinction_probabilities().x

    print(f'average growth rate {l} with stable proportions {v}')
    #print(f'Monte Carlo simulated extinction rates: {simulated_extinction_rates}')
    #print(f'How close (holistic) to a zero: {distance_from_zero}')
    print(f'True extinction probability: {true_extinction_rates}')
    plot.plot_extinction_monte_carlo_graph(simulated_extinction_rates, true_extinction_rates, 
                                   simulations_per_type=simulations_per_type)