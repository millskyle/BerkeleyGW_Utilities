#!/usr/bin/python

f = open("POSCAR",'r')
o = open("input",'w')


infile = []
for line in f:
   infile.append(line.split())



class System(object):
   def comment(self,data):
      self.comment = data.join()



s = System()
s.comment = infile[0]


print s.comment







