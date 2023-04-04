import constants
import models

if __name__ == '__main__':
    #params = constants.TOY_DYING_PARAMS
    params = constants.TOY_LIVING_PARAMS
    model = models.NoncollaborativeStochastic(params)
    det_model = models.NoncollaborativeDeterministic(params)
    l, v = det_model.get_stable_state()
    
    simulated_extinction_rates = model.sample_extinction(1000)
    distance_from_zero = model.get_extinction_function()(simulated_extinction_rates)
    true_extinction_rates = model.find_extinction_probabilities(simulated_extinction_rates).x

    print(f'average growth rate {l} with stable proportions {v}')
    print(f'Monte Carlo simulated extinction rates: {simulated_extinction_rates}')
    print(f'How close (holistic) to a zero: {distance_from_zero}')
    print(f'True extinction probability: {true_extinction_rates}')