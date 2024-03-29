import re

# # Opening file and extracting and assigning accessions to a list
# file_path = r"D:\Uni\Data Science\P4B_exercises\exercises and examples\regular_expressions\exercises\accessions.txt"
# file = open(file_path, 'r')
# file_1 = file.read()
#
# accessions = re.findall(r"'([^']+)'", file_1)
#
# contains_5 = []
# contain_d_or_e = []
# contain_d_and_e = []
# contain_d_x_e = []
# contain_d_e = []
# start_x_or_y = []
# start_x_or_y_end_e = []
# contain_num3row = []
# end_d_or = []
#
# for ac_no in accessions:
#     if matches:= re.search("[5]", ac_no):
#         contains_5.append(ac_no)
#     if matches:= re.search("[de]", ac_no):
#         contain_d_or_e.append(ac_no)
#     if matches := re.search("de", ac_no):
#         contain_d_and_e.append(ac_no)
#     if matches := re.search("d.e", ac_no):
#         contain_d_x_e.append(ac_no)
#     # if matches := re.search("de|ed", ac_no):
#     #     contain_d_e.append(ac_no)
#     if matches := re.search(r"d.*e|e.*d", ac_no):
#         contain_d_e.append(ac_no)
#     if matches := re.search("^[xy]", ac_no):
#         start_x_or_y.append(ac_no)
#     if matches := re.search(r"^[xy].*e$", ac_no):
#         start_x_or_y_end_e.append(ac_no)
#     if matches := re.search("\d{3,}", ac_no):
#         contain_num3row.append(ac_no)
#     if matches := re.search("d[arp]$", ac_no):
#         end_d_or.append(ac_no)
#
#
# print(f"Contains 5: {contains_5}")
# print(f"Contains the letter d or e: {contain_d_or_e}")
# print(f"Contains the letters d and e in that order: {contain_d_and_e}")
# print(f"Contains the letters d and e in that order with a single letter between: {contain_d_x_e}")
# print(f"Contains both d and e in any order: {contain_d_e}")
# print(f"Start with x or y: {start_x_or_y}")
# print(f"Start with x or y and end with e: {start_x_or_y_end_e}")
# print(f"Contain 3 or more numbers in a row: {contain_num3row}")
# print(f"End with d followed by either a, r or p: {end_d_or}")

f = open('dna.txt','r')
dna = f.readline()
print(dna)
fragments = []

pattern = r"(\w*A[AGTC]T)(AAT\w*)"
pattern_1 = r"(\w*GC[AG][AT])(TG\w*)"


while len(dna) != 0:

    frag = re.split(f"{pattern}[AGTC]*{pattern_1}", dna)
    print(frag)
    if len(frag) < 2:
        fragments.append(frag[0])
        break
    dna = frag[1]
    fragments.append(frag[2])

print(fragments[::-1])
