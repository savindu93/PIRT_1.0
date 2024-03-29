f = open("seq.txt", 'r')

# for line in f:
#     print(line.strip())

# string = f.read()

species = {}

# for str in f:
#     list = str.split()
#     species[list[0]] = [list[1] + " " + list[2], list[3]]
#
# print(species)

species_1 = {}
for str in f:
    list = str.split()
    species_1[list[0]] = {'Scientific name': list[1] + " " + list[2], 'Anatomy':list[3]}


print(species_1)

print(species_1['Rice']['Scientific name'])






# list = string.split()
# print(string)
# print(list)





