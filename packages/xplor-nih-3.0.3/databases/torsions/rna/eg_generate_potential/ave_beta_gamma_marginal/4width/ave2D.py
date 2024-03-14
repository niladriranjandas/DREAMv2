
import glob


# First energy file.
infile = open(glob.glob('../../beta_gamma_delta/4width_10grid/beta_gamma_marginal/*_energy.dat')[0], 'r')
contents = infile.read()
infile.close()

# Read name, atom, info, axis, and energy lines.
atomlines = []
infolines1 = []
axislines = []
energy1 = []
for line in contents.splitlines():
    if line.startswith('name'): name1 = line.split()[1]
    if line.startswith('atom'): atomlines.append(line)
    if line.startswith('info'): infolines1.append(line)
    if line.startswith('axis'): axislines.append(line)
    if line.startswith('energy'):
        words = line.split()
        for word in words[1:]:
            energy1.append(float(word))        
        


# Second energy file.
infile = open(glob.glob('../../alpha_beta_gamma/4width_10grid/beta_gamma_marginal/*_energy.dat')[0], 'r')
contents = infile.read()
infile.close()

# Read name, info, and energy lines.
infolines2 = []
energy2 = []
for line in contents.splitlines():
    if line.startswith('name'): name2 = line.split()[1]
    if line.startswith('info'): infolines2.append(line)
    if line.startswith('energy'):
        words = line.split()
        for word in words[1:]:
            energy2.append(float(word))


# Average energy values.
energy = []
for (e1, e2) in zip(energy1, energy2):
    energy.append((e1+e2)/2.0)



# Output file.
filename = '_'.join(name1.split('_')[:3]) + '_energy.dat'
outfile = open(filename, 'w')

# Write name and atom lines.
outfile.write('name %s\n\n' % '_'.join(name1.split('_')[:3]))
for line in atomlines:
    outfile.write(line+'\n')

# Write info lines
outfile.write('\n')
outfile.write('info Average of %s and %s\n' % (name1, name2))
outfile.write('info Info on %s\n' % name1)
for line in infolines1:
    outfile.write(line+'\n')
outfile.write('info Info on %s\n' % name2)
for line in infolines2:
    outfile.write(line+'\n')
outfile.write('\n')

# Write axis lines.
for line in axislines:
    outfile.write(line+'\n')
outfile.write('\n')

# Write energy lines.
linesize = 10
count = 0
for value in energy:
    if count == 0:
        outfile.write('energy %.3f ' % value)
        count += 1
    elif count > 0 and count < linesize-1:
        outfile.write('%.3f ' % value)
        count +=1
    else:
        outfile.write('%.3f \n' % value)
        count = 0
outfile.write('\n\n\n')

outfile.close()








