class animal():

    def __init__(self, nombre, edad, ruido):
        self.nombre = nombre
        self.edad = edad
        self.ruido = ruido

    def haz_ruido(self):
            print(self.ruido)


ave = animal("piolin", 2, "pio pio")
ave.haz_ruido()
print(ave.__dict__)


class perro(animal):
    raza = ''


class gato(animal):
    usa_arenero = True


perro = perro("Rocky", 4, "guau guau")
perro.raza = "rottweiler"
print(perro.__dict__)

gato = gato("michi", 3, "miau miau")
print(gato.__dict__)
