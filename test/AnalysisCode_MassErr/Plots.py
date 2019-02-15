#!/usr/bin/env python

import os
import glob
import sys

PathsToPlots = [['Results/Type0/', 'Type0'], ['Results/Type2/', 'Type2'], ['../UsefulScripts/StabilityCheck/pictures/', 'StabilityPlots'], ['../UsefulScripts/DeDxStudy74X/pictures_FAST/', 'DeDxStudyPlots'], ['../UsefulScripts/DeDxStudy74X/systematics/', 'DeDxSystematics/']]
FilesLists = []

def produceIndex_html (listOfPlots):
   f = open('index.html', 'w')
   f.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">\n')
   f.write('<html>\n')
   f.write('<head>\n')
   f.write('<title>Analysis Plots  \'s Images</title>\n')
   f.write('</head>\n')
   f.write('<BODY TEXT="#000066" BGCOLOR="#CCCCFF"  LINK="#FFFF00" VLINK="#FF00FF" >\n')
   f.write('<P>\n')
   f.write('<TABLE WIDTH="100%">\n')
   f.write('<TR>\n')
   f.write('<TH ALIGN="center" width="40%" BGCOLOR="#9999FF"><FONT \n')
   f.write('COLOR="#FFFFFF" SIZE=+1> Analysis Plots --- Type 0 --- </FONT></TH>\n')
   f.write('</TABLE>\n')
   f.write('</p>\n\n\n')
   f.write(' <center>\n')
   f.write('<TABLE border=0 cellpadding=5 cellspacing=2>\n')
   for i in range(0, len(listOfPlots), 3):
      f.write('<TR> <TD align=center> <a href="%s"><img src="%s"hspace=5 vspace=5 border=0 style="width: 80%%" ALT="%s"></a>\n' % (os.path.basename(listOfPlots[i]), os.path.basename(listOfPlots[i]), os.path.basename(listOfPlots[i])))
      f.write('  <br> %s </TD>\n' % os.path.basename(listOfPlots[i]))
      if i+1 < len(listOfPlots):
         f.write('  <TD align=center> <a href="%s"><img src="%s"hspace=5 vspace=5 border=0 style="width: 80%%" ALT="%s"></a>\n' % (os.path.basename(listOfPlots[i+1]), os.path.basename(listOfPlots[i+1]), os.path.basename(listOfPlots[i+1])))
         f.write('  <br> %s </TD>\n' % os.path.basename(listOfPlots[i+1]))
      if i+2 < len(listOfPlots):
         f.write('  <TD align=center> <a href="%s"><img src="%s"hspace=5 vspace=5 border=0 style="width: 80%%" ALT="%s"></a>\n' % (os.path.basename(listOfPlots[i+2]), os.path.basename(listOfPlots[i+2]), os.path.basename(listOfPlots[i+2])))
         f.write('  <br> %s </TD> </TR>\n' % os.path.basename(listOfPlots[i+2]))

   f.write('</TABLE>\n')
   f.write(' </center>\n')
   f.write('<p>\n')
   f.write('<hr width="100%">\n')
   f.write('<table border=0 cellspacing=0 cellpadding=0 width="100%">\n')
   f.write('<tr><td><em>\n')
   f.write('Created: %s\n\n' % os.popen('date').read())
   f.write('</em></td><td align=right><em>\n\n')
   f.write('</table>\n')
   f.write('</body>\n')
   f.write('</html>\n')
   f.close()

def getListOfPlots ():
   for directory in PathsToPlots:
      FilesLists.append(glob.glob('%s*.png' % directory[0]))

if sys.argv[1] == '1':
   getListOfPlots()
   CWD = os.getcwd()
   for i in range (0, len(PathsToPlots)):
      os.system('mkdir -p PlotsSummary/%s' % PathsToPlots[i][1])
      for file in FilesLists[i]:
         os.system('cp %s PlotsSummary/%s' % (file, PathsToPlots[i][1]))
      os.chdir('PlotsSummary/%s' % PathsToPlots[i][1])
      produceIndex_html (FilesLists[i])
      os.chdir(CWD)

