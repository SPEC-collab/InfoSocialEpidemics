# Information/Social Epidemiology
#
# SPEC collaborative
#
# Key elements to remember
#
# 1. All specifics are coded into mobility
# 2. Two grids are needed
# 3. Find what the earmark means for our agents (Graham scheduler)
#
# 2021

from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
from datacollection import DataCollector
from scipy.stats import poisson, bernoulli
from grahamscheduler import GrahamActivation
from enum import Enum
import numpy as np
import random
import sys


# Strip down to SIRD (susceptible, infected, recovered, deceased)

class Stage(Enum):
    SUSCEPTIBLE = 1
    INFECTED = 2
    RECOVERED = 3
    DECEASED = 4

# For mobility:
#
# ISOLATED: limited mobility (i.e. low probability of moving to another place in the same grid)
# -> mobility_prob = dist(central, dispersion) ~ very sharp Poisson distribution
#
# LOCAL_ONLY: average mobility (i.e. 0.5 probability of moving within the same grid)
# -> mobility_prob = P(X) ~ 0
#
# BETWEEN_GRIDS: probability of moving within the local grid + prob going to the other grid = 1
# -> switch_grid_prob = 1 - prob(mobility_prob)
# -> local_dwell_time
# -> intercity_dwell_time

class MobilityType(Enum):
    ISOLATED = 1
    LOCAL_ONLY = 2
    BETWEEN_GRIDS = 3

class ISEpiAgent(Agent):
    """ An agent representing a potential covid case"""
    
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.stage = Stage.SUSCEPTIBLE
       
        # These are fixed values associated with properties of individuals
        self.incubation_time = poisson.rvs(model.avg_incubation)
        self.dwelling_time = poisson.rvs(model.avg_dwell)
        self.recovery_time = poisson.rvs(model.avg_recovery)
        self.prob_contagion = self.model.prob_contagion
       
        # Mortality in vulnerable population appears to be around day 2-3
        self.mortality_value = self.model.mortality_rate/(self.model.dwell_15_day*self.recovery_time)
        # Severity appears to appear after day 5
        self.severity_value = self.model.prob_severe/(self.model.dwell_15_day*self.recovery_time)
        self.curr_dwelling = 0
        self.curr_incubation = 0
        self.curr_recovery = 0
        
    def alive(self):
        print(f'{self.unique_id} {self.age_group} {self.sex_group} is alive')

    def is_contagious(self):
        return (self.stage == Stage.EXPOSED) or (self.stage == Stage.ASYMPTOMATIC)

    # In this function, we count effective interactants
    def interactants(self):
        count = 0

        if (self.stage != Stage.DECEASED) and (self.stage != Stage.RECOVERED):
            for agent in self.model.grid.get_cell_list_contents([self.pos]):
                if agent.unique_id != self.unique_id:
                    if not(agent.isolated) or self.isolated_but_inefficient:
                        count = count + 1

        return count

    def step(self):
        # Using the model, determine if a susceptible individual becomes infected due to
        # being elsewhere and returning to the community
        if self.stage == Stage.SUSCEPTIBLE:
            if bernoulli.rvs(self.model.rate_inbound):
                self.stage = Stage.INFECTED

        if self.stage == Stage.SUSCEPTIBLE:
            # Important: infected people drive the spread, not
            # the number of healthy ones
            cellmates = self.model.grid.get_cell_list_contents([self.pos])
            infected_contact = False

            # Isolated people should only be contagious if they do not follow proper
            # shelter-at-home measures
            for c in cellmates:
                if c.is_contagious():
                    c.add_contact_trace(self)
                    if self.isolated and bernoulli.rvs(1 - self.model.prob_isolation_effective):
                        self.isolated_but_inefficient = True
                        infected_contact = True
                        break
                    else:
                        infected_contact = True
                        break        
            
            if infected_contact:
                if bernoulli.rvs(self.prob_contagion):
                    self.stage = Stage.INFECTED

            # Second opportunity to get infected: residual droplets in places
            # TODO

            if not(self.isolated):
                self.move()
        elif self.stage == Stage.INFECTED:
            if self.curr_incubation + self.curr_recovery < self.incubation_time + self.recovery_time:
                self.curr_recovery = self.curr_recovery + 1

                if bernoulli.rvs(self.severity_value):
                    self.stage = Stage.SEVERE
            else:
                self.stage = Stage.RECOVERED
        elif self.stage == Stage.RECOVERED:
            cellmates = self.model.grid.get_cell_list_contents([self.pos])
            
            # A recovered agent can now move freely within the grid again
            self.curr_recovery = 0
            self.move()
        elif self.stage == Stage.DECEASED:
            pass
        else:
            # If we are here, there is a problem 
            sys.exit("Unknown stage: aborting.")

        self.astep = self.astep + 1

    def move(self):
        # If dwelling has not been exhausted, do not move
        if self.curr_dwelling > 0:
            self.curr_dwelling = self.curr_dwelling - 1

        # If dwelling has been exhausted, move and replenish the dwell
        else:
            possible_steps = self.model.grid.get_neighborhood(
                self.pos,
                moore=True,
                include_center=False
            )
            new_position = self.random.choice(possible_steps)

            self.model.grid.move_agent(self, new_position)
            self.curr_dwelling = poisson.rvs(self.model.avg_dwell)

def compute_susceptible(model):
    return count_type(model, Stage.SUSCEPTIBLE)/model.num_agents

def compute_infected(model):
    return count_type(model, Stage.INFECTED)/model.num_agents

def compute_recovered(model):
    return count_type(model, Stage.RECOVERED)/model.num_agents

def compute_deceased(model):
    return count_type(model, Stage.DECEASED)/model.num_agents

def count_type(model, stage):
    count = 0

    for agent in model.schedule.agents:
        if agent.stage == stage:
            count = count + 1

    return count

def compute_isolated(model):
    count = 0

    for agent in model.schedule.agents:
        if agent.isolated:
            count = count + 1

    return count/model.num_agents

def compute_employed(model):
    count = 0

    for agent in model.schedule.agents:
        if agent.employed:
            count = count + 1

    return count/model.num_agents

def compute_unemployed(model):
    count = 0

    for agent in model.schedule.agents:
        if not(agent.employed):
            count = count + 1

    return count/model.num_agents

def compute_contacts(model):
    count = 0

    for agent in model.schedule.agents:
        count = count + agent.interactants()

    return count/len(model.schedule.agents)

def compute_stepno(model):
    return model.stepno

def compute_eff_reprod_number(model):
    prob_contagion = model.prob_contagion
    infected = 0.0
    infected_time = 0.0

    for agent in model.schedule.agents:
        if agent.stage == Stage.INFECTED:
            infected += 1
            infected_time = infected_time + agent.incubation_time

        else:
            continue

    avg_contacts = compute_contacts(model)
    return prob_contagion * avg_contacts * (infected_time/infected)

def compute_num_agents(model):
    return model.num_agents

class ISEpiModel(Model):
    """
    A model to understand the effect of social and information aspects as effects that
    change the course of events for the usual SIRD model due to changes in the underlying
    mean field theory.
    """
    dwell_15_day = 96
    
    def __init__(self, epidemiology, gridworld, network, social, scheduler, dummy=0):
        self.running = True
        self.schedule = GrahamActivation(self) if (scheduler == "graham") else RandomActivation(self)
        self.stepno = 0
    
        # Dwell times per day and average number of them per agent
        self.dwell_15_day = 96
    
        # Setup of the grids
        self.num_agents_a = gridworld["city_a"]["num_agents"]
        self.num_agents_b = gridworld["city_b"]["num_agents"]
        
        self.commuters_a = gridworld["city_a"]["commuters"]
        self.commuters_b = gridworld["city_b"]["commuters"]
        
        self.grid_a = MultiGrid(gridworld["city_a"]["width"], gridworld["city_a"]["height"], True)
        self.grid_b = MultiGrid(gridworld["city_b"]["width"], gridworld["city_b"]["height"], True)
        
        # Epidemiology-related parameters
        self.avg_incubation = int(round(epidemiology["avg_incubation_time"] * self.dwell_15_day))
        self.prob_contagion = epidemiology["prob_contagion"]
        self.avg_dwell = epidemiology["avg_dwell"]
        self.avg_recovery = epidemiology["avg_recovery_time"] * self.dwell_15_day
        self.num_init = int(self.num_agents * epidemiology["prop_initial_infected"])
        self.mortality_rate = epidemiology["mortality_rate"]

        # Setup city A
        
        self.i = 0

        for ag in self.age_distribution:
            for sg in self.sex_distribution:
                r = self.age_distribution[ag]*self.sex_distribution[sg]
                num_agents = int(round(self.num_agents*r))
                mort = self.age_mortality[ag]*self.sex_mortality[sg]
                for k in range(num_agents):
                    a = ISEpiAgent(self.i, ag, sg, mort, self)
                    self.schedule.add(a)
                    x = self.random.randrange(self.grid.width)
                    y = self.random.randrange(self.grid.height)
                    self.grid.place_agent(a, (x,y))
                    self.i = self.i + 1
                    
        # Setup city B
        
        self.datacollector = DataCollector(
            model_reporters = {
                "Step": compute_stepno,
                "NA": compute_num_agents_a,
                "NB": compute_num_agents_b,
                "SusceptibleA": compute_susceptible_a,
                "InfectedA": compute_infected_a,
                "RecoveredA": compute_recovered_a,
                "DeceasedA": compute_deceased_a,
                "IsolatedA": compute_isolated_a,
                "LocalOnlyA": compute_local_only_a,
                "BetweenCitiesA": compute_between_cities,
                "RtA": compute_eff_reprod_number_a,
                "SusceptibleB": compute_susceptible_b,
                "InfectedB": compute_infected_b,
                "RecoveredB": compute_recovered_b,
                "DeceasedB": compute_deceased_b,
                "IsolatedB": compute_isolated_b,
                "LocalOnlyB": compute_local_only_b,
                "BetweenCitiesB": compute_between_cities_b,
                "RtB": compute_eff_reprod_number_b
            }
        )

        # Final step: infect an initial proportion of random agents  
        for a in self.schedule.agents:
            if self.num_init < 0:
                break
            else:
                a.stage = Stage.EXPOSED
                num_init = num_init - 1
   
    def step(self):
        self.datacollector.collect(self)
        
        if self.stepno % self.dwell_15_day == 0:
            print(f'Simulating day {self.stepno // self.dwell_15_day}')
        
        self.schedule.step()
        self.stepno = self.stepno + 1
