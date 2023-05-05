import numpy as np
from scipy.special import comb
from scipy.linalg import expm
from scipy.optimize import root, fixed_point
from scipy.integrate import odeint




class Model:
  def __init__(self, params):
    self.params = params

  def generate_timepoint_data(self, n_initial, duration : float, timestep_count : int, stochastic=True, normalized=False):
    timestep = duration / timestep_count
    timepoints = [timestep * i for i in range(timestep_count + 1)]
    data = [n_initial]
    current = n_initial
    for _ in range(timestep_count):
      if stochastic:
        current = self.run(current, timestep)
      else:
        current = self.run_det(current, timestep)
      data.append(current)
    if normalized:
      _, growth = self.calculate_quasistable_distribution()
      scaled_data = []
      for t, row in zip(timepoints, data):
        factor = np.exp(t * growth)    
        normalized_row = [datapoint / factor for datapoint in row]
        scaled_data.append(normalized_row)
      data = scaled_data    
    return np.array(data), timepoints

  def run(self, n_initial, t, max_steps = 100000):
    return n_initial

  def run_det(self, n_initial, t):
    return n_initial

  def sample_extinction(self, attempts_per_type : int):
    M = self.params['M']
    extinction_rates = []
    t = 100
    for i in range(M + 1):
      print(f'running {attempts_per_type} simulations starting with single cell of type {i}')
      extinct_count = 0
      for _ in range(attempts_per_type):
        try: 
          n_initial = [0] * (M + 1)
          n_initial[i] = 1
          outcome = self.run(n_initial, t, max_steps=10000)
          if sum(outcome) == 0:
            extinct_count += 1
        except RateTooLargeException:
          continue
      extinction_rates.append(extinct_count / attempts_per_type)
    return extinction_rates 

  def calculate_quasistable_distribution(self):
    return 0, 0


class LinearModel(Model):
  """To instantiate a Linear Model, you will need to implement the following properties:
  - all rate methods: _r_b(), _r_d(), _r_mu(), _r_um()
  - the number of sites: self.M
  """

  def __init__(self, params):
    super().__init__(params)
    self._initialize_events()
    self._initialize_transtion_matrix()

  def run(self, n_initial, t, max_steps = 100000):
    n = n_initial[:]
    step_count = 0
    while True:
      step_count += 1
      # don't run forever! 
      if step_count > max_steps:
        raise RateTooLargeException(f'simulation exceeded {max_steps} steps!')
      # calculate rate of incoming events
      total_rate = self._get_total_rate(n)
      # if total rate is 0, there will be no more events
      if total_rate == 0:
        break
      # sample time until next event
      waiting_time = - np.log(np.random.random()) / total_rate
      # update time remaining after next event
      t -= waiting_time
      # if time remaining is negative, next event came too late
      if t < 0:
        break
      # sample which event is coming next
      event_rng = total_rate * np.random.random()
      for event in self.events:
        event_rate = event["rate"] * n[event["cell"]]
        event_rng -= event_rate
        if event_rng < 0:
          found_event = event
          break
      # update cell counts by what the event does
      self._implement_event(found_event, n)
    # simulation completed: either time ran out, or the population died out.
    return n

  def run_det(self, n_initial, t):
    return np.array(n_initial) @ expm(self.flow_matrix * t)
    
  def calculate_extinction_rates(self):
    M = self.params['M']
    def function_to_solve(x):
      return np.subtract(self._update_extinction_guess(x), x)
    guess = [min(self._r_d(i / M) / self._r_b(i / M), 1) for i in range(M + 1)]
    return root(function_to_solve, guess)

  def calculate_quasistable_distribution(self):
    M = self.params['M']
    def function_to_solve(x):
      return self._update_quasistable_guess(x)
    unnormalized_guess = self.run_det([1 / (M + 1) for i in range(M + 1)], 10)
    s = sum(unnormalized_guess)
    guess = [g / s for g in unnormalized_guess]
    dist = fixed_point(function_to_solve, guess, maxiter = 2000)
    growth = np.log(sum(self.run_det(dist, 1)))
    return dist, growth

  def calculate_limit_extinction_rates(self, point_count = 1000):
    
    def extinction_derivative(y, x):
      if methylation_evolution(x) == 0:
        if y == 1:
          return 0
        return (self._r_d(x) - self._r_b(x)) / (self._r_um(x) + self._r_mu(x))
      return ((self._r_b(x) * y - self._r_d(x)) * (y - 1)) / (self._r_mu(x) * x - self._r_um(x) * (1 - x))
  
    def methylation_evolution(i):
      return self._r_mu(i) * i - self._r_um(i) * (1 - i)
    
    initial_x = root(methylation_evolution, 0.5).x[0]
    initial_y = min(self._r_d(initial_x) / self._r_b(initial_x), 1)
    print(f'initial extinction: {initial_y} at {initial_x}')
    all_times = [i / point_count for i in range(point_count + 1)]
    low_times = [time for time in all_times if time < initial_x] + [initial_x]
    high_times = [initial_x] + [time for time in all_times if time > initial_x]



    

    # get probabilities below the "stable" x
    low_res = odeint(extinction_derivative, initial_y, list(reversed(low_times)))

    # get probabilities above the "stable" x  
    high_res = odeint(extinction_derivative, initial_y, high_times)

    # combine probabilities into one list
    res = list(reversed(low_res[:,0][1:])) + list(high_res[:,0][1:])

    return res

  def _update_extinction_guess(self, initial_guess):
    new_guess = []
    M = self.params['M']
    p = self.params['p']
    for i in range(M + 1):
      numerator = 0
      denominator = 0
      methylation_level = i / M
      # demethylation
      if i != 0:
        rate = self._r_mu(methylation_level) * i
        numerator += rate * initial_guess[i - 1]
        denominator += rate
      # methylation
      if i != M:
        rate = self._r_um(methylation_level) * (M - i)
        numerator += rate * initial_guess[i + 1]
        denominator += rate 
      # birth
      rate = self._r_b(methylation_level)
      numerator += rate * (sum([comb(i, j) * (p ** j) * ((1 - p) ** (i - j)) * initial_guess[j] for j in range(i + 1)]) ** 2)
      denominator += rate
      # death
      rate = self._r_d(methylation_level)
      numerator += rate
      denominator += rate
      new_guess.append(numerator / denominator)
    return new_guess

  def _update_quasistable_guess(self, initial_guess):
    unscaled = self.run_det(initial_guess, 1)
    new_guess = unscaled / sum(unscaled)
    return new_guess
  
  def _initialize_events(self):
    self.events = []
    M = self.params['M']
    for i in range(M + 1):
      # include the events for a cell with i sites methylated
      self.events.append({'type': 'birth', 'rate': self._r_b(i / M), 'cell': i})
      self.events.append({'type': 'death', 'rate': self._r_d(i / M), 'cell': i})
      if i != M:
        self.events.append({'type': 'methylation', 'rate': (M - i) * self._r_um(i / M), 'cell': i})
      if i != 0:
        self.events.append({'type': 'demethylation', 'rate': i * self._r_mu(i / M), 'cell': i})

  def _initialize_transtion_matrix(self):
    M = self.params['M']
    flow_matrix = []
    for i in range(M + 1):
      flow_row = []
      methylation_level = i / M
      outflow = 0 # flow out of type i cell
      for j in range(M + 1):
        # calculate the flow from i into j
        inflow = 0 # clow into type j cell
        # if demethylation causes flow
        if i - 1 == j:
          increment = self._r_mu(methylation_level) * i
          inflow += increment
          outflow += increment
        # if methylation causes flow
        if i + 1 == j:
          increment = self._r_um(methylation_level) * (M - i)
          inflow += increment
          outflow += increment
        # if death causes flow
        if i == j:
          increment = self._r_d(methylation_level)
          outflow += increment
        # if birth causes inflow
        if i >= j:
          p = self.params['p']
          increment = 2 * (p ** j) * (1 - p) ** (i - j) * comb(i, j) * self._r_b(methylation_level)
          inflow += increment
        # if birth causes outflow
        if i == j:
          increment = self._r_b(methylation_level)
          outflow += increment
        flow_row.append(inflow)
      # to calculate total flow matrix, subtract off the outward flow from type i cell
      flow_row[i] = flow_row[i] - outflow
      flow_matrix.append(flow_row)
    self.flow_matrix = np.array(flow_matrix)


  def _r_b(self, i):
    return self.params['b_0'] * (1 - i) + i * (self.params['b_M'])

  def _r_d(self, i):
    return self.params['d_0'] * (1 - i) + i * (self.params['d_M'])
  
  def _r_um(self, i):
    return self.params['r_um']
  
  def _r_mu(self, i):
    return self.params['r_mu']
  
  def _get_total_rate(self, n):
    subtotal_rate = 0
    for event in self.events:
      cell_type = event['cell']
      subtotal_rate += event['rate'] * n[cell_type]
    return subtotal_rate
  
  def _implement_event(self, event, n):
    event_type = event["type"]
    cell = event["cell"]
    n[cell] = n[cell] - 1
    if event_type == 'birth':
      p = self.params['p']
      for _ in range(2):
        new_cell = np.random.binomial(cell, p)
        n[new_cell] = n[new_cell] + 1 
    elif event_type == 'methylation':
      n[cell + 1] = n[cell + 1] + 1 
    elif event_type == 'demethylation':
      n[cell - 1] = n[cell - 1] + 1 
    elif event_type == 'death':
      pass
    else: 
      print('unrecognized event type')

class ThresholdModel(LinearModel):

  def _r_b(self, i):
    c = self.params['c']
    if i > c: 
      return self.params['b_M']
    return self.params['b_0']
  
  def _r_d(self, i):
    c = self.params['c']
    if i > c: 
      return self.params['d_M']
    return self.params['d_0']

class CollaborativeModel(LinearModel):
  def __init__(self, params, degree=1):
    self.degree = degree
    super().__init__(params)


  def _r_mu(self, i):
    return self.params['r_mu'] + self.params['r_mu_u'] * (1 - i) ** self.degree
  
  def _r_um(self, i):
    return self.params['r_um'] + self.params['r_um_m'] * i ** self.degree

class SuperCollaborativeModel(LinearModel):
  
  def _r_mu(self, i):
    return self.params['r_mu'] + self.params['r_mu_u'] * (1 - i) ** 2
  
  def _r_um(self, i):
    return self.params['r_um'] + self.params['r_um_m'] * i ** 2





class RateTooLargeException(Exception):
  """Raised when the total rate is too large to be worthwhile for simulations to continue"""
