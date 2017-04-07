#!/usr/bin/env python

"""

.. module:: fortran2py-Module
    :synopsis: This module attempts to automatize the fortran to python code translation\
    (specific to boltzmann code) 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import re,sys

fortranfile = 'gravitinoBRs.f'

ffile = open(fortranfile,'r')
fcode = ffile.read()
ffile.close()

Dic = {"dlog" : "log", "dexp" : "exp", "dsqrt" : "sqrt", "zeta3" : "Zeta3"
       , "pi" : "Pi", "mp" : "MP", "dble" : "float", ".gt." : " > ", ".lt." : " < ",
       ".eq." : " == ", ".and." : " and ", ".or." : " or ", ".ne." : " != "}

#Fix all comments:
fcode = fcode.replace("\nC","\n#")
fcode = fcode.replace("\nc","\n#")
fcode = fcode.replace("\n!","\n#")
fcode = fcode.replace("!","#")
#Transform everything to lower case
fcode = fcode.lower()
#Put together same line
patt = re.compile("(\n     [\$\&\,])")
matches =  patt.findall(fcode)
for match in matches: fcode = fcode.replace(match,"")
#Remove trailing blanks spaces
fcode = fcode.replace("\n      ","\n")
#Replace dictionary:
for key in Dic: fcode = fcode.replace(key,Dic[key])
#Remove double precision notation
patt = re.compile("(\d\d*\.*\d*d[-+\d]\d*)")
matches =  patt.findall(fcode)
for match in matches:
    replacement = match.replace("d","*10.**(")
    replacement += ")"
    if "*10.**(0)" in replacement:
        if "." in match: replacement = replacement.replace("*10.**(0)","")
        else: replacement = replacement.replace("*10.**(0)",".")
    fcode = fcode.replace(match,replacement)

#Split lines
fcode = fcode.split("\n")
#Remove GOTO trailing numbers and add comment
goto_labels = []
for il,line in enumerate(fcode):
    patt = re.compile("goto[\d]*")
    cline = line.replace(" ","")
    for match in patt.findall(cline): goto_labels.append(match.replace("goto",""))
    goto_labels = list(set(goto_labels))
    for label in goto_labels:
        if line[:len(label)] == label:
#             fcode[il] = "#The line below contained the GOTO label "+ label + "\n" +line.replace(label,"") 
            fcode[il] = line.replace(label,"")
#Remove right spaces:
for il,line in enumerate(fcode): fcode[il] = line.rstrip()
#Remove variable declaration and store variable names:
vars_dec = ["double precision","integer","external","logical","real","common","save"
            ,"character","double complex", "double coMPlex"]
var_names = []
for il,line in enumerate(fcode):
    if not line:
        fcode[il] = None
        continue
    for v in vars_dec:        
        if line[:len(v)] == v:            
            if v != "common" and v != "save":
                line = line.replace(v,"")
                line = line[line.find(" "):].lstrip().replace(" ","")
            else:       
                line = line[line.rfind("/")+1:].lstrip().replace(" ","")
            if not "," in line and line:
                var_names.append(line)
            else:
                var_names.append(line[:line.find(",")])
                var_names.append(line[line.rfind(",")+1:])
                matches = re.compile(",[\w]*,").findall(line)
                matches += re.compile(",[\w]*\([\d,]*\),").findall(line)
                for im,match in enumerate(matches): matches[im] = match.lstrip(',').rstrip(',')
                var_names += matches
                var_names = list(set(var_names))
            var_names += line.split(",")
            fcode[il] = None            
while fcode.count(None) > 0: fcode.remove(None)
#Switch fortran array notation to python:
var_arrays =  []
for var in var_names:
    if "(" in var and ")" in var: var_arrays.append(var[:var.index("(")])  
for il,line in enumerate(fcode):
    for var in var_arrays:
        while True:
            if not var+"(" in line: break
            i1 = line.find(var+"(")
            oldvar = line[i1:]
            i2 = oldvar.find(")")            
            oldvar = oldvar[:i2+1]
            newvar = oldvar
            while "," in newvar: newvar = newvar.replace(",","][")
            newvar = newvar.replace("(","[").replace(")","]")
            line = line.replace(oldvar,newvar)
            
        fcode[il] = line

#Fix identation
ident = 0
for il,line in enumerate(fcode):
    patt = re.compile("^do[\w]*=[\d],[\d]|^do[\w]*=[\d],[\d],[\d]|^if\(.*\)then|dowhile")    
    cline = line.replace(" ","")
    if cline == "enddo" or cline == "endif": ident -= 4
    if cline[0:4] == "else" or cline[0:6] == "elseif": newline = " "*(ident-4) + line.lstrip()
    else: newline = " "*ident + line.lstrip()
    if patt.search(cline): ident += 4        
    fcode[il] = newline
    
#Convert if and do blocks:
for il,line in enumerate(fcode):
    cline = line.replace(" ","")
    if cline == "enddo" or cline == "endif":
        fcode[il] = "\n"
        continue  
    patt1 = re.compile("^if\(.*\)then")
    patt2 = re.compile("do[\w]*=[\d],[\d]|^do[\w]*=[\d],[\d],[\d]")
    cline = line.replace(" ","")
    if patt1.search(cline):
        newline = line.replace("if(","if (").lstrip("(")
        newline = newline.rstrip("then").rstrip(")")
        newline += " :"
    elif patt2.search(cline):
        newline = line.replace('do','for')
        newline = newline.replace("=", " in range(")
        end = newline[newline.find(" in range(")+10:].replace(" ","") + "):"
        newline = newline[:newline.find(" in range(")+10] + end
    else: newline = line
    fcode[il] = newline
         

pyfile = fortranfile.replace(".f",".py").replace(".F",".py")
outfile = open(pyfile,"w")    
for line in fcode: outfile.write(line+"\n")
outfile.close()
print 'Converted',fortranfile,'to',pyfile
