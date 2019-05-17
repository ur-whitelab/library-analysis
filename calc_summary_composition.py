import sys
import operator

if(len(sys.argv) != 2):
    print('Usage: calc_summary_composition.py [target_AA_strings_file]. Must be formatted correctly.')

fname = sys.argv[1]

with open(fname, 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    lines[i] = lines[i].replace('\n','')

print(lines)

summary_dict = {}
for line in lines:
    if line not in summary_dict.keys():
        summary_dict[line] = 1
    else:
        summary_dict[line] += 1

sorted_summary = sorted(summary_dict.items(), key = operator.itemgetter(1), reverse=True)

print(sorted_summary[:100])

with open('{}_SUMMARY.txt'.format(fname.split('.')[0]), 'w+') as f:
    for item in sorted_summary:
        f.write('{}, {}\n'.format(item[0], item[1]))
