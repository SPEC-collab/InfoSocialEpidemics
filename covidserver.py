# Information/Social Epidemiology
#
# SPEC collaborative
#
# 2021

# A simple tunable model for COVID-19 response
from mesa.visualization.modules import CanvasGrid
from mesa.visualization.modules import ChartModule
from mesa.visualization.ModularVisualization import ModularServer
from mesa.visualization.UserParam import UserSettableParameter

from isepi import ISEpiModel
from isepi import Stage
from isepi import MobilityType

def agent_portrayal(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "Layer": 0,
                 "r": 0.5}

    if agent.stage == Stage.SUSCEPTIBLE:
        portrayal["Color"] = "blue"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.EXPOSED:
        portrayal["Color"] = "red"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.ASYMPTOMATIC:
        portrayal["Color"] = "brown"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.SYMPDETECTED:
        portrayal["Color"] = "yellow"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.ASYMPDETECTED:
        portrayal["Color"] = "cyan"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.SEVERE:
        portrayal["Color"] = "magenta"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.RECOVERED:
        portrayal["Color"] = "green"
        portrayal["Layer"] = 0
    elif agent.stage == Stage.DECEASED:
        portrayal["Color"] = "black"
        portrayal["Layer"] = 0
    elif agent.locked:
        portrayal["Color"] = "gray"
        portrayal["Layer"] = 0
    else:
        portrayal["Color"] = "gray"
        portrayal["Layer"] = 0

    return portrayal

grid = CanvasGrid(agent_portrayal, 50, 50, 400, 400)

chart = ChartModule([{"Label": "Susceptible",
                      "Color": "Blue"},
                      {"Label": "Infected",
                      "Color": "Red"},
                      {"Label": "Recovered",
                      "Color": "Green"},
                      {"Label": "Deceased",
                      "Color": "Black"}
                      ],
                    data_collector_name='datacollector')

chart_epidemiology = ChartModule([
                      {"Label": "Rt",
                      "Color": "Blue"
                      },
                      ],
                    data_collector_name='datacollector'
)

chart_number = ChartModule([
                      {"Label": "N",
                      "Color": "Blue"
                      },
                      ],
                    data_collector_name='datacollector'
)

model_params = {
    "num_agents": 240,
    "width": 50,
    "height": 50,
    "prop_initial_infected": UserSettableParameter("slider", "Proportion of initially infected", 0.001, 0.0, 0.1, 0.001),
    "avg_incubation_time": UserSettableParameter("slider", "Average incubation time", 5, 2, 24, 1),
    "avg_recovery_time": UserSettableParameter("slider", "Average recovery time", 15, 15, 30, 1),
    "mortality_rate": UserSettableParameter("slider", "Mortality rate", 0.13, 0.0, 0.20, 0.01),
    "prob_contagion": UserSettableParameter("slider", "Probability of contagion", 0.03, 0.0, 0.1, 0.005),
}

server = ModularServer(ISEpiModel,
                       [grid,chart,chart_epidemiology, chart_number],
                       "Information- and social-based epidemic model",
                       model_params
                       )

server.port = 8521 # The default
server.launch()