import numpy as np
from scipy.special import comb
from scipy.linalg import expm, eig
from scipy.optimize import root
from scipy.integrate import odeint
class Model:

  def generate_timepoint_data(self, n_initial, duration : float, timestep_count : int):
    timestep = duration / timestep_count
    timepoints = [timestep * i for i in range(timestep_count + 1)]
    data = [n_initial]
    current = n_initial
    for _ in range(timestep_count):
      current = self.run(current, timestep)
      data.append(current)
    return np.array(data), timepoints

  def run(self, n_initial, t):
    return n_initial

class NoncollaborativeStochastic(Model):
  def __init__(self, params):
    self.M = params[0]
    self.b = params[1]
    self.d_0 = params[2]
    self.d_M = params[3]
    self.r_um = params[4]
    self.r_mu = params[5]
    self.p = params[6]
    
    self.events = []
    for i in range(self.M + 1):     # births
      self.events.append({"type": 0, "rate": self.b, "start": i})
    for i in range(self.M + 1):     # deaths
      self.events.append({"type": 1, "rate": self.d_0 + (self.d_M - self.d_0) * i / self.M, "start": i})
    for i in range(self.M):         # methylations
      self.events.append({"type": 2, "rate": self.r_um * (self.M - i), "start": i})
    for i in range(1, self.M + 1):  # demethylations
      self.events.append({"type": 3, "rate": self.r_mu * i, "start": i})
  
  def run(self, n_initial, t, max_steps = 100000):
    n = [_ for _ in n_initial]
    step_count = 0
    while True:
      step_count += 1
      # get total rate
      total_rate = 0
      for event in self.events:
        total_rate += self._get_rate(n, event)

      if step_count > max_steps:
        raise RateTooLargeException('Simulation exceeded maximum steps!')


      # get waiting time
      if total_rate == 0:
        break
      waiting_time = - np.log(np.random.random()) / total_rate
      t -= waiting_time
      if t <= 0:
        break

      # get event
      event_RNG = np.random.random() * total_rate
      for event in self.events:
        event_RNG -= self._get_rate(n, event)
        if event_RNG <= 0:
          break
        
      # implement event
      self._implement_event(event, n)

    return n

  
  def sample_extinction(self, attempts_per_type : int):
    extinction_rates = []
    t = 100
    for i in range(self.M + 1):
      print(f'running {attempts_per_type} simulations starting with single cell of type {i}')
      extinct_count = 0
      for _ in range(attempts_per_type):
        try: 
          n_initial = [0] * (self.M + 1)
          n_initial[i] = 1
          outcome = self.run(n_initial, t, max_steps=10000)
          if sum(outcome) == 0:
            extinct_count += 1
        except RateTooLargeException:
          continue
      extinction_rates.append(extinct_count / attempts_per_type)
    return extinction_rates 

  def find_extinction_probabilities(self):
    # we calculate the probabilities if all r_mu and r_um are 0
    guess = [min(self._r_d(i) / self._r_b(i), 1) for i in range(self.M + 1)]
    f = self._get_extinction_function()
    return root(f, guess)

  def find_limit_extinction_probabilities(self, point_count = 1000):
    func = self._get_limit_extinction_derivative()
    initial_x = self.r_um / (self.r_mu + self.r_um)
    initial_y = self._get_limit_extinction_initial_probability()
    all_times = [i / point_count for i in range(point_count + 1)]
    low_times = [time for time in all_times if time < initial_x] + [initial_x]
    high_times = [initial_x] + [time for time in all_times if time > initial_x]

    # get probabilities below the "stable" x
    low_res = odeint(func, initial_y, list(reversed(low_times)))

    # get probabilities above the "stable" x  
    high_res = odeint(func, initial_y, high_times)

    # combine probabilities into one list
    res = list(reversed(low_res[:,0][1:])) + list(high_res[:,0][1:])

    return res




  def _get_extinction_function(self):
    return lambda x : np.subtract(x, self._get_extinction_by_first_step(x))

  def _get_limit_extinction_derivative(self):
    def derivative(y, x):
      if x == self.r_um / (self.r_mu + self.r_um):
        if y == 1:
          return 0
        return ((self._r_d(self.M * x) - self._r_b(self.M * x))
            / (self._r_d(self.M * x) * self._r_b(self.M * x) * (self.r_um + self.r_mu)))
      return (self._r_b(x * self.M) * y - self._r_d(x * self.M)) * (1 - y) \
        / (- x * self.r_mu + (1 - x) * self.r_um) 
    return derivative

  def _get_limit_extinction_initial_probability(self):
    x = (self.r_um) / (self.r_mu + self.r_um)
    return min(self._r_d(self.M * x) / self._r_b(self.M * x), 1)

  def _r_d(self, i):
    return (self.d_0 * (self.M - i) + self.d_M * i) / self.M

  def _r_m(self, i):
    return (self.M - i) * self.r_um
  
  def _r_u(self, i):
    return i * self.r_mu

  def _r_b(self, i):
    return self.b

  def _f(self, i : int, event_type : str):
    if event_type == 'b':
      num = self._r_b(i)
    elif event_type == 'd':
      num = self._r_d(i)
    elif event_type == 'u':
      num = self._r_u(i)
    elif event_type == 'm':
      num = self._r_m(i)
    else:
      raise TypeError('event_type should be one of b, d, u, or m')
    
    denom = self._r_b(i) + self._r_d(i) + self._r_u(i) + self._r_m(i)
    
    return num / denom

  def _get_rate(self, n, event):
      return n[event["start"]] * event["rate"]


  def _get_extinction_by_first_step(self, extinction_rates):
    new_extinction_rates = []
    for i in range(self.M + 1):
      p = 0
      if i != 0:
        p += self._f(i, 'u') * extinction_rates[i - 1]
      if i != self.M:
        p += self._f(i, 'm') * extinction_rates[i + 1]
      p += self._f(i, 'b') * (extinction_rates[i] ** 2)
      p += self._f(i, 'd')
      new_extinction_rates.append(p)
    return new_extinction_rates
  

  def _implement_event(self, event, n):
    event_type = event["type"]
    event_start = event["start"]
    if event_type == 0:   # birth
      self._implement_birth(event_start, n)
    elif event_type == 1: # death
      self._implement_death(event_start, n)
    elif event_type == 2: # methylation
      self._implement_methylation(event_start, n)
    elif event_type == 3: # demethylation
      self._implement_demethylation(event_start, n)

  def _implement_birth(self, n_m, n):
    n[n_m] = n[n_m] - 1
    for _ in range(2):
      n_m_new = 0
      for __ in range(n_m):
        if np.random.random() <= self.p:
          n_m_new += 1
      n[n_m_new] = n[n_m_new] + 1
  
  def _implement_death(self, n_m, n):
    n[n_m] = n[n_m] - 1

  def _implement_methylation(self, n_m, n):
    n[n_m] = n[n_m] - 1
    n[n_m + 1] = n[n_m + 1] + 1

  def _implement_demethylation(self, n_m, n):
    n[n_m] = n[n_m] - 1
    n[n_m - 1] = n[n_m - 1] + 1

class NoncollaborativeDeterministic(Model):
  def __init__(self, params):
    self.M = params[0]
    self.b = params[1]
    self.d_0 = params[2]
    self.d_M = params[3]
    self.r_um = params[4]
    self.r_mu = params[5]
    self.p = params[6]

    # make evolution matrix
    T = []
    for i in range(self.M + 1):
      T.append([self._get_flow(i, j) for j in range(self.M + 1)])    
    self.T = np.array(T)
  
  def run(self, n_initial, t):
    n = np.array(n_initial) @ expm(self.T * t)
    return n

  def get_stable_state(self):
    w, vl = eig(self.T, left=True, right=False)
    
    # find max eigenvalue
    for i, k in enumerate(w):
      if i == 0:
        candidate = k
        index = i
      if k > candidate:
        candidate = k
        index = i
      
    state = vl[:, index]
    return np.real(candidate), state / np.sum(state)
    
  # calculates the flow from i -> j
  def _get_flow(self, i, j):
    subtotal = 0    
    if i == j - 1:  # inward flow due to methylation
      subtotal += (self.M - i) * self.r_um
    if i == j + 1:  # inward flow due to demethylation
      subtotal += i * self.r_mu
    if i >= j:      # inward flow due to birth
      subtotal += 2 * self.b * comb(i, j) * (self.p) ** j * (1 - self.p) ** (i - j)
    if i == j:      # outward flows
      subtotal -= (self.d_0 + i / self.M * (self.d_M - self.d_0))  # death
      subtotal -= i * self.r_mu              # demethylation
      subtotal -= (self.M - i) * self.r_um   # methylation
      subtotal -= (self.b)                   # birth
    return subtotal


class NoncollaborativeStochasticFull(Model):
  """This is a refactored, full version of the noncollaborative model that takes into account
  varying birthrates and a better way (dictionary) to specify parameters"""

  def __init__(self, params):
    self.b_0 = params['b_0']
    self.b_M = params['b_M']
    self.d_0 = params['d_0']
    self.d_M = params['d_M']
    self.r_um = params['r_um']
    self.r_mu = params['r_mu']
    self.p = params['p']
    self.M = params['M']

    self.events = []

    for i in range(self.M + 1):
      self.events.append({'type': 'birth', 'rate': self._r_b(i), 'cell': i})
      self.events.append({'type': 'death', 'rate': self._r_d(i), 'cell': i})
      self.events.append({'type': 'methylation', 'rate': self._r_um(i), 'cell': i})
      self.events.append({'type': 'demethylation', 'rate': self._r_mu(i), 'cell': i})

  def run(self, n_initial, t, max_steps = 100000):
    n = n_initial[:]
    while True:

      total_rate = self._get_total_rate(n)
      # if total rate is 0, there will be no more events
      if total_rate == 0:
        break
      waiting_time = - np.log(np.random.random()) / total_rate
      t -= waiting_time
      # if next time of event is too late, no more events
      if t < 0:
        break

      # find which event is next
      event_rng = total_rate * np.random.random()
      for event in self.events:
        event_rate = event["rate"] * n[event["cell"]]
        event_rng -= event_rate
        if event_rng < 0:
          found_event = event
          break
      # update cell counts by what the event does
      self._implement_event(found_event, n)
    return n

  def _r_b(self, i):
    return self.b_0 + (self.b_M - self.b_0) * i / self.M

  def _r_d(self, i):
    return self.d_0 + (self.d_M - self.d_0) * i / self.M
  
  def _r_um(self, i):
    return (self.M - i) * self.r_um
  
  def _r_mu(self, i):
    return self.M * self.r_mu
  
  def _get_total_rate(self, n):
    subtotal_rate = 0
    for event in self.events:
      subtotal_rate += event['rate'] * n[event['cell']]
    return subtotal_rate
  
  def _implement_event(self, event, n):
    event_type = event["type"]
    cell = event["cell"]
    if event_type == 'birth':
      n[cell] = n[cell] + 1
    elif event_type == 'death':
      n[cell] = n[cell] - 1
    elif event_type == 'methylation':
      n[cell] = n[cell] - 1
      n[cell + 1] = n[cell + 1] + 1 
    elif event_type == 'demethylation':
      n[cell] = n[cell] - 1
      n[cell - 1] = n[cell - 1] + 1 
    else: 
      print('unrecognized event type')

class RateTooLargeException(Exception):
  """Raised when the total rate is too large to be worthwhile for simulations to continue"""
