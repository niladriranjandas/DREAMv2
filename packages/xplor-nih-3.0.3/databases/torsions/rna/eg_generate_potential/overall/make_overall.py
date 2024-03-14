
import glob

outfile = open('rna09_v0.dat', 'w')


infile = open(glob.glob('../delta-1_epsilon-1_zeta-1/4width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../epsilon-1_zeta-1_alpha/4width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../zeta-1_alpha_beta/4width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../alpha_beta_gamma/4width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../beta_gamma_delta/4width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../delta_chi/2width_10grid/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../ave_epsilon-1_zeta-1_marginal/4width/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../ave_zeta-1_alpha_marginal/4width/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../ave_alpha_beta_marginal/4width/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../ave_beta_gamma_marginal/4width/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')

infile = open(glob.glob('../ave_delta/2width_4width/*dat')[0], 'r')
contents = infile.read()
infile.close()

outfile.write(contents+'\n\n')



