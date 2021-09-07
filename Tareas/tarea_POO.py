"""
NAME
       tarea_POO.py
VERSION
        0.0.1
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This is the example of a class to create characters in a video game,
        use inheritance to create villains and heroes.
CATEGORY
        Video game
GITHUB REPOSITORY
        https://github.com/rod13-afk/python_class/blob/master/Tareas/tarea_POO.py
"""


# This is a class to create the characters of a video game
class character():

    def __init__(self, name, gender):
        self.name = name
        self.gender = gender
        self.attack = 10
        self.defense = 10
        self.max_speed = 8
        self.speed = 0
        self.health = 13
        self.stamina = 9

    def move_up(self):
        self.speed = self.max_speed - 1

    def move_down(self):
        self.speed = self.max_speed + 1

    def move_right(self):
        self.speed = self.max_speed

    def move_left(self):
        self.speed = self.max_speed

    def attack(self):
        self.stamina = self.stamina - 1

    def defend(self):
        self.health = self.health - 1


class villain(character):

    def __init__(self, name, gender):
        super().__init__(name, gender)
        self.rob_a_bank = True

    def move(self):
        pass


class hero(character):
    def __init__(self, name, gender):
        super().__init__(name, gender)
        self.defeat_the_villain = True

    def move_up(self):
        pass


villain = villain("Chucho", "Hombre")
print(villain.__dict__)

hero = hero("Mateo", "Hombre")
print(hero.__dict__)
