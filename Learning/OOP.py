class Dog:

    attr1 = "Mammal"

    def __init__(self, name):

        self.name = name


rodger = Dog("Rodger")
print(f"Rodger is a {Dog.attr1}")

class person:

    def __init__(self, name, id):

        self.name = name
        self.id = id
        self.__ent = "Person"

    def details(self):
        print(f"Person's name is {self.name}\n"
              f"Id No.: {self.id}\n"
              f"Entity: {self.__ent}\n")


class employee(person):

    def __init__(self, name, id, salary, post):

        super().__init__(name, id)
        self.salary = salary
        self.post = post

    def details(self):
        print(f"Employee's name is {self.name}\n"
              f"Id No.: {self.id}\n"
              f"Salary: {self.salary}\n"
              f"Post: {self.post}\n")


p1 = person("Nimal", 1098)
p2 = employee("Rahul", 2098, 200000, "SE")

p1.details()
p2.details()
print(f"{p2.name} {p2.id}")



