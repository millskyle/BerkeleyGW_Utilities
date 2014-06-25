#!/usr/bin/python
exclude=[]
import sys

class read_POSCAR(object):
   global exclude
   x=0; y=1; z=2; lvec=[]
   def __init__(self,filename):
      x=0;y=1;z=2
      self.infile = []
      for line in open(filename,'r'):
         self.infile.append(line.split())    #fill the infile list with the lines from the file
      self.comment = " ".join(self.infile[0])  #first line is a comment
      self.lc = float(" ".join(self.infile[1])) #second line is the lattice constant (lc)
      for i in [2,3,4]:
         self.lvec.append([float(j) for j in self.infile[i]])  #Next three lines are lattice vectors
      self.atoms = self.infile[5] #sixth line is a list atom types
      self.natoms = 0 #define an attribute to hold the total number of atoms
      self.count_atoms = { 'atomic symbol':'number of type' } #define a dict to hold counts and atomic symbols
      self.atom_type_list=[]  # blank list to hold atom type list
      for i,atom in enumerate(self.atoms): #For each *type* of atom,
         if not atom in exclude: #if we're not excluding this type of atom
            self.count_atoms[self.atoms[i]] = int(self.infile[6][i]) #add the count to the dict
            for _ in range(int(self.infile[6][i])): #add $count instances of the atomic symbol
               self.atom_type_list.append(self.atoms[i]) #to the symbol list
            self.natoms+=int(self.infile[6][i]) #add the count of this type to the total atom count
      self.ntyp = len(set(self.atom_type_list)) #get the number of different atom types
      if self.infile[7][0][0].upper()=="S": #if selective dynamics
         self.infile.pop(7)  #remove the line that says "selective dynamics" (keeps the line numbers in sync)
      if self.infile[7][0][0].upper() =="D": #if direct coordinates
         self.coord_scaling=[ self.lvec[0][x], self.lvec[1][y], self.lvec[2][z] ]   #then all coordinates 
                                                                                    #must be multiplied by 
                                                                                    #the lattice vectors
      elif self.infile[7][0][0].upper() =="C": #if cartesian, then we're going to multply by 1.0.
         self.coord_scaling=[ 1.,1.,1. ]
      else: #if something else, then I don't know what to do so just leave them as-is
         self.coord_scaling=[ 1.,1.,1. ]

      self.position_data = [] #define list to hold all position data
      for i in range(8,8+self.natoms): #start at "line 8" and loop over natoms lines
         self.position_data.append([
            self.atom_type_list[i-8], #the atomic symbol
            float(self.infile[i][x]) * self.coord_scaling[x], #the x-coordinate, scaled
            float(self.infile[i][y]) * self.coord_scaling[y], #the y-coordinate, scaled
            float(self.infile[i][z]) * self.coord_scaling[z] ])  #the z-coordinate, scaled

def out(data):
   global o
   #Simple function to write data to file and print it to stdout.
   o.write(data+"\n")
   print data

def main():
   global exclude
   global o
   if len(sys.argv)<3:
      print """
ERROR!
      USAGE:
            python POSCAR2PWscf.py inputPOSCAR outputFile "<atoms to exclude>"

       (where <atoms to exclude> is a space separated list of atom types to exclude from the PWSCF input file,
       ie:  "C O" would exclude all carbon and oxygen atoms present in the POSCAR from the output)
"""
      ThrowError()

   if len(sys.argv)>3:
      exclude=sys.argv[3].split()
   else:
      exclude=[]

   x=0;y=1;z=2
   o = open(sys.argv[2],'w') #open the output file for writing
   s = read_POSCAR('POSCAR') #initialize the POSCAR object (ie: read it and do the processing)

   #start printing the file out.  Options are hard-coded.
   out("&control")
   out("""   prefix = '!PREFIX!'
   calculation = 'scf'
   restart_mode = 'from_scratch'
   wf_collect = .true.
   tstress = .true.
   tprnfor = .true.
   outdir = './'
   wfcdir = './'
   pseudo_dir = './'""")
   out("/")
   out("&system")
   out("   celldm(1) = {0}".format(s.lc))
   out("""   ibrav = 0
   nbnd = 60
   ecutwfc = 160.0""")
   out("   nat = {0}".format(s.natoms))
   out("   ntype = {0}".format(s.ntyp))
   out("/")
   out("&electrons")
   out("""   electron_maxstep = 100
   conv_thr = 1.0d-10
   mixing_mode = 'plain'
   mixing_beta = 0.7
   mixing_ndim = 8
   diagonalization = 'david'
   diago_david_ndim = 4
   diago_full_acc = .true.""")
   out("/")
   out("CELL PARAMETERS")
   for i in s.lvec:
      out("\t{0}\t{1}\t{2}".format(i[0],i[1],i[2]))

   out("ATOMIC SPECIES")
   for i in s.atoms:
      if not i in exclude:
         out("\t{0}\t{1}\t{2}.UPF ".format(i,atomic_masses[i],i))
   out("ATOMIC_POSITIONS angstrom")
   for i in s.position_data:
      if not i[0] in exclude:
         out("\t{0}\t{1}\t{2}\t{3}".format(i[0],i[1],i[2],i[3]))
   out("...")
   out("K_POINTS gamma")
   out("")


#An atomic mass list...verify these before running.
atomic_masses = {
"H":1.0079,
"He":4.0026,
"Li":6.941,
"Be":9.0122,
"B":10.811,
"C":12.0107,
"N":14.0067,
"O":15.9994,
"F":18.9984,
"Ne":20.1797,
"Na":22.9897,
"Mg":24.305,
"Al":26.9815,
"Si":28.0855,
"P":30.9738,
"S":32.065,
"Cl":35.453,
"Ar":39.948,
"K":39.0983,
"Ca":40.078,
"Sc":44.9559,
"Ti":47.867,
"V":50.9415,
"Cr":51.9961,
"Mn":54.938,
"Fe":55.845,
"Co":58.9332,
"Ni":58.6934,
"Cu":63.546,
"Zn":65.39,
"Ga":69.723,
"Ge":72.64,
"As":74.9216,
"Se":78.96,
"Br":79.904,
"Kr":83.8,
"Rb":85.4678,
"Sr":87.62,
"Y":88.9059,
"Zr":91.224,
"Nb":92.9064,
"Mo":95.94,
"Tc":98,
"Ru":101.07,
"Rh":102.9055,
"Pd":106.42,
"Ag":107.8682,
"Cd":112.411,
"In":114.818,
"Sn":118.71,
"Sb":121.76,
"Te":127.6,
"I":126.9045,
"Xe":131.293,
"Cs":132.9055,
"Ba":137.327,
"La":138.9055,
"Ce":140.116,
"Pr":140.9077,
"Nd":144.24,
"Pm":145,
"Sm":150.36,
"Eu":151.964,
"Gd":157.25,
"Tb":158.9253,
"Dy":162.5,
"Ho":164.9303,
"Er":167.259,
"Tm":168.9342,
"Yb":173.04,
"Lu":174.967,
"Hf":178.49,
"Ta":180.9479,
"W":183.84,
"Re":186.207,
"Os":190.23,
"Ir":192.217,
"Pt":195.078,
"Au":196.9665,
"Hg":200.59,
"Tl":204.3833,
"Pb":207.2,
"Bi":208.9804,
"Po":209,
"At":210,
"Rn":222,
"Fr":223,
"Ra":226,
"Ac":227,
"Th":232.0381,
"Pa":231.0359,
"U":238.0289,
"Np":237,
"Pu":244,
"Am":243,
"Cm":247,
"Bk":247,
"Cf":251,
"Es":252,
"Fm":257,
"Md":258,
"No":259,
"Lr":262,
"Rf":261,
"Db":262,
"Sg":266,
"Bh":264,
"Hs":277,
"Mt":268
}



main()
