class Species:


    __international_code = ''
    n_objects = 0

    def __init__(self, generic_name, species_name, TID = "NA", num_protein = "NA", PID = "NA"):

        self.generic_name = generic_name
        self.sname = species_name
        self.TID = TID
        self.num_protein = num_protein
        self.PID = PID

        Species.n_objects += 1


    def getFullName(self):
        return f"{self.generic_name} {self.sname}"

    @classmethod
    def updateCode(cls, code):
        cls.__international_code = code

    @classmethod
    def getIntCode(cls):
        return f"{cls.__international_code}"

    @staticmethod
    def welcomeMessage():
        print("Welcome to The Species Database")

class Subspecies(Species):

    def __init__(self, generic_name, species_name, TID = "NA", num_protein = "NA", PID = "NA",
                 subname = "NA", g_type = "NA"):

        super().__init__(generic_name, species_name, TID, num_protein, PID)

        self.subname = subname
        self.g_type = g_type

    def getFullName(self):
        return f"{self.generic_name} {self.sname} subsp.{self.subname}"





if __name__ == '__main__':

    p1 = Species('Zea', 'mays', '4577', 330923, '01')
    p2 = Species('Oryza', 'sativa')

    p3 = Subspecies('Oryza', 'sativa', '39947', 330923, '03', 'japonica', 'circular')
    # p4 = Subspecies('Oryza', 'sativa', '39948', 'indica', 'long')

    print(f"{p3.generic_name} {p3.sname} {p3.subname}\n"
          f"Protein ID: {p3.PID}")


    print(p3.getFullName())

    p3.updateCode("ICNafp")
    print(p3.getIntCode())



    # print(f"{p4.generic_name} {p4.sname} {p4.subname}\n"
    #       f"Protein ID: {p4.PID}")



