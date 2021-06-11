# SPEC collaborative
#
# Graham Scheduler code based on MESA
#
# 2021

from mesa.time import BaseScheduler
from mesa.model import Model


class GrahamActivation(BaseScheduler):
    """A scheduler devised by Jeff Graham the simulates concurrent, random activation of 
    pockets of agents. Agents have 'earmarks', identifiers that are associated with local
    ordering of events within one single timestep.
    """

    def __init__(self, model: Model, partitions: int):
        """Create an empty Staged Activation schedule.
        Args:
            model: Model object associated with the schedule.
            partitions: number of partitions indicating how agents 
            can be divided into ordered classes.
        """
        super().__init__(model)
        self.stages = {}
        self.partitions = partitions


    def step(self) -> None:
        """Executes the step of all agents according to their earmarks, randomly
        for each earmark.
        """

        for p in range(0, self.partitions):
            self.stages[p] = []

        for agent in self._agents.values(): 
            self.stages[agent.earmark].append(agent.unique_id)
           
        for earmark in self.stages.keys():
            self.model.random.shuffle(self.stages[earmark])

            for unique_id in self.stages[earmark]:
                self._agents[unique_id].step()

        self.steps += 1
        self.time += 1