from abc import ABC, abstractmethod


class AbstractCommand(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def execute(self):
        pass
