import constants
import models

if __name__ == '__main__':
    #params = constants.TOY_DYING_PARAMS
    params = constants.TOY_LIVING_PARAMS
    model = models.NoncollaborativeStochastic(params)
    det_model = models.NoncollaborativeDeterministic(params)
    l, v = det_model.get_stable_state()
    print(l, v)
    extinction_rates = model.get_extinction(1000)
    print(extinction_rates)