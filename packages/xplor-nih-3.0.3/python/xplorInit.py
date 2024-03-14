
"""python interface intialization script

this should not normally be modified by users.

this executes in module ``xplor''
"""

import __main__

import os, sys

#set locale and encoding. This can't be allowed to change all over the place.
import locale
locale.setlocale(locale.LC_ALL,("C","utf8"))

#
# set up argv, scriptName, progName
#
scriptName=""
progName=""
try:
  progName=os.environ["PROGNAME"]
except KeyError:
  pass

#reset thse Python env vars, so that launched programs don't inherit them
import os
for k in os.environ.keys():
  if k.startswith('PYTHON'):
    del os.environ[k]
    pass
  pass

#generalize for DYLD_LIBRARY_PATH, remove only paths which include
#XPLORNIHROOT
import os
XPLOR_DIR=os.environ['XPLOR_DIR']
LD_LIBPATH_VAR = os.environ['LD_LIBPATH_VAR']
import xplorSimulation
xplorSimulation.orig_LD_LIBRARY_PATH=os.environ[LD_LIBPATH_VAR]
os.environ["orig_"+LD_LIBPATH_VAR]=os.environ[LD_LIBPATH_VAR]

#
# 2020/05/06 - LD_LIBRARY_PATH should have different values depending on
# whether it is used for binaries in this distribution or execed programs.
# An example of one is the web-browser lauched by the xplorJupyter notebook
# command.
#

paths=xplorSimulation.orig_LD_LIBRARY_PATH.split(':')
newPath=":".join([path for path in paths if not path.startswith(XPLOR_DIR)])
os.environ[LD_LIBPATH_VAR]=newPath

#
##
## rework importing to also reset the library load path
##
#import builtins
#savedImport=builtins.__import__
#import os
#def xplorImport(name,lglobals=None, llocals=None, fromlist=(), level=0):
#  if not lglobals: lglobals=globals()
#  if not llocals:  llocals=locals()
#  savedVal = os.environ[LD_LIBPATH_VAR]
#  os.environ[LD_LIBPATH_VAR] = os.environ["orig_"+LD_LIBPATH_VAR]
#  print("importing...",name," using ",os.environ[LD_LIBPATH_VAR])
#  retVal = savedImport(name,lglobals, llocals, fromlist, level)
#  os.environ[LD_LIBPATH_VAR] = savedVal 
#  return retVal
#  pass
#builtins.__import__=xplorImport
#

#
# scriptName is either the first or the last argument - this will be
# executed by the Python interpreter if it's present.
#
cnt=1
from sys import argv
newargv=argv[1:]

if len(argv)>1:
  arg1 = argv[1]
  argN = argv[-1]
  if arg1=='-m':
    newargv = argv[3:]
    scriptName = argv[3]
    progName=argv[2]
  elif not arg1.startswith('-'):
    # first argument - rewrite argv to be everything else
    newargv = argv[2:]
    scriptName = arg1
    progName=arg1
  elif not argN.startswith('-'):
    # last argument - remove from argv
    newargv = argv[:-1]
    scriptName = argN
    pass
  pass
sys.argv = [progName] + newargv


import sys
if sys.version_info[0] >= 3:
  def execfile(filename, globals=None, locals=None):
    if globals is None:
      globals = sys._getframe(1).f_globals
      pass
    if locals is None:
      locals = sys._getframe(1).f_locals
      pass
    with open(filename, "r") as fh:
      code = compile(fh.read(), filename, 'exec')
      exec(code, globals, locals)
      pass
    return
  pass
else:
  execfile = execfile
  pass
      

#
# Usually stderr and stdout should be line-buffered. When there is no tty,
# these care fully block-buffered by default, and this leads to much
# confusion, and possibly lost output.
#
lineBuffering=True
import os ; varName = "XPLOR_NOLINEBUFFER"
if varName in os.environ and os.environ[varName]=="1":
  lineBuffering=False
  pass
  
import io
import sys
for out in "stdout stderr".split():
  tmp=eval("sys.%s.detach()" % out )
  if hasattr(tmp,"detach"): tmp = tmp.detach()
  exec("sys.%s = io.TextIOWrapper(tmp,line_buffering=lineBuffering)" % out)
  pass
#
#tmp=sys.stdout.detach() #same for stdout - very confusing if stderr comes early
#if hasattr(tmp,"detach"):
#  tmp = tmp.detach()
#sys.stdout = io.TextIOWrapper(tmp,
#                              line_buffering=True) 
#

def parseArguments(options=[],
                   cmdline="",
                   description="[no description provided]",
                   usageString="",
                   allowExtraOptions=False):
  """process extra Xplor-NIH command-line arguments.

  Mistyped options will be detected.

  Additional options can be specified using the options argument.
  Options start with a single or double dash. Option can have one or more
  arguments, the number of which is specified by a number separated by a
  colon from the option name. For option arguments the literal character
  sequence %q% is converted to a single quote. Thus, it is not possible to
  pass this sequence as an option argument.

  If allowExtraOptions is True, extra options will be allowed, and not treated
  as errors.   

  Options present on the command-line but which are not specified in the
  options argument will cause an error.

  The return value is (optlist,args), where
    optlist is [(option, arg1, arg2..), ...]
    and args are the remaining command-line arguments after options have
    been parsed.

  example:
    (optlist,args) = parseArguments(["help-script","outfile:1"])

  The cmdline, description and usageString arguments are used to generate output
  when the following options are specified:

      --help-cmdline  - prints the program name (as deduced from argv[0]),
                        followed by the cmdline argument.
      --help-description - print the description string with newline characters
                           replaced with spaces.
      --help-script      - print a usage string build from the cmdline,
                           description and usageString arguments.
      --output-rst       - the usage string printed using --help-script 
                           normally has its reStructedText markup stripped
                           with the <m rst2txt> facility. If this option is
                           placed first, the markup is preserved.
      --pdf              - in a similar fashion to the --output-rst option,
                           if this option is placed before --help-script, the
                           reStructuredText help string is converted to PDF
                           and displayed in a separate window. This option
                           has the side-effect of setting the module-global
                           variable outputPDF to True (otherwise it is False).

  For each of these options after printing the string, exit() is called.
  """

  if type(options)==type("string"):
    options = (options,)

  import getopt # only for its exception

  # first parse the options list
  optionList={}
  optionsString=""
  for option in options:
    numArgs=0
    if ':' in option:
      (option,numArgs) = option.split(':')
      pass
    numArgs = int(numArgs)
    optionList[option] = numArgs

    optionsString += "  -%s" % option
    if numArgs:
      delim='='
      for cnt in range(numArgs):
        optionsString += "%s%s" % (delim,"value")
        delim=','
        pass
      pass
    optionsString += '\n'

    pass


  from sys import argv
  import os.path
  name = os.path.basename(argv[0])
  cmdline = name + " " + cmdline
  cnt = 1
  retOptions=[]
  retArgs=[]
  noMoreOpts=False
  outputRST=False
  global outputPDF
  outputPDF=False
  doHelpScript=False
  try:
    while cnt < len(argv):
      arg = argv[cnt]
      cnt += 1
      if arg.startswith('-') and not noMoreOpts:
        if arg=='--':
          noMoreOpts=True
          continue
        option = arg.lstrip('-')
        option_vals = option.split('=')
        option = option_vals[0]
        if len(option_vals)>1:
          vals = option_vals[1].split(',')
        else:
          vals=None
          pass
        if option=="help-description":
          print(description.replace('\n',' ').strip().replace('  ',' '))
          exit(0)
        if option=="help-cmdline":
          print(cmdline.replace('\n',' '))
          exit(0)
        if option=="output-rst":
          outputRST=True
          continue
        if option=="pdf":
          outputPDF=True
          continue
        if option=="help-script" and not option in optionList.keys():
          doHelpScript=True
          continue
        if not option in optionList.keys():
          raise getopt.GetoptError("unrecognized option: %s"%arg)
        num=int(optionList[option])
          
        entry=[option]
        if vals:
          if num!=len(vals):
            raise getopt.GetoptError("wrong number of values for option: " +
                                     option)
          entry += vals
          pass
        else:
          for cnt2 in range(num):
            if cnt >= len(argv):
              raise getopt.GetoptError("not enough values for option: " +
                                       option);
            entry.append( argv[cnt].replace('%q%',"'") )
            cnt += 1
            pass
          pass
        retOptions.append(entry)
        pass
      else:
        retArgs.append(arg)
        pass
      pass
    pass
  except getopt.GetoptError as exc:
    if not allowExtraOptions:
        global scriptName
        if progName==scriptName: scriptName=""
        writeConsole("\n%s\n\n" % exc.msg)
        writeConsole("Usage: %s\n" % cmdline)
        writeConsole(" where options is zero or more of\n")
        writeConsole(optionsString)
        writeConsole("\nor an Xplor-NIH command option.\n")
        writeConsole(\
          "To list the Xplor-NIH options specify --help on the command line.\n")
        import sys
        sys.exit(2)
        pass
    pass
  if doHelpScript:
    import sys
    inputString=''
    if not sys.stdin.isatty():
      inputString=sys.stdin.read()
    
    from docutils.core import publish_string
    rstString = name + '\n'
    rstString += "-"*len(name) + '\n'
    rstString += '\n' + description + '\n'
    rstString += "Usage: " + cmdline + '\n'
    rstString += usageString
    if inputString:
      rstString += r'''
.. raw:: latex

  \clearpage

'''
      pass
    rstString += inputString
    import os
    dir=os.path.join(os.environ['XPLOR_HELPDIR'],'helperScripts')

    from utils import printReStructuredText
    printReStructuredText(rstString,
                          outputType="rst" if outputRST else
                          "pdf" if outputPDF else "txt",
                          cachePrefix=name,
                          cacheDir=dir,
                          cacheTime=0)
    exit(0)
  return (retOptions,retArgs)
  

def requireVersion(vStr,
                   raiseException=True):
  """Returns True if the current Xplor-NIH version is less than
  that specified in vStr. If raiseException=True, an exception is
  raised.
  """
  import xplor
  vNums=vStr.split('.')

  import types
  for i in range(len(vNums)):
    reqNum=int(vNums[i])
    vNum = version_info[i] if i<len(version_info) else 0
    if type(vNum)!=type(1): vNum=0
    if reqNum<vNum: break
    if reqNum>vNum:
      if raiseException:
        raise exception("Xplor-NIH version is too old. Please upgrade.")
      else:
        return False
    pass
  return True

def rcDir():
  """return the path of the Xplor-NIH configuration directory.

  If the directory does not exist, create it
  """
  import os

  ret = os.path.join(os.environ['HOME'],'.xplornih')
  if os.path.exists(ret) and not os.path.isdir(ret):
    raise exception("%s is not a directory")

  return ret

def rcPyDir():
  """return the path of the Xplor-NIH configuration directory.

  If the directory does not exist, create it
  """
  import os

  ret = os.path.join(rcDir(),'python')
  if os.path.exists(ret) and not os.path.isdir(ret):
    raise exception("%s is not a directory")

  return ret

def pyConfigFile(name):
  """return the full path of a .xplornih configuration file whose filename
  is name
  """

  import os
  return os.path.join(rcPyDir(),name)


from os import environ as env
def help():
  return open(env['XPLOR_HELPDIR'] + '/nih-py-xplor').read()
def fastCommand(c,r=()):
   """call <m xplorWrap>.XplorWrap.fastCommand
   """
   return wrap.fastCommand(c,r)


#
# fix up system path - remove /lib/python2.x components of strings
#
tmp=[]
import sys, os
for e in sys.path:
  e=e.replace('/dist/lib/python2.5','/dist')
  e=e.replace('/dist/lib/python2.4','/dist')
  e=e.replace('/dist/lib/python2.3','/dist')

  # remove cwd and nonexistent entries
  try:
    if not os.path.samefile(os.path.curdir,e):
      tmp.append(e)
      pass
  except OSError:
    pass
  

sys.path = tmp

sys.path.append( rcPyDir() )
sys.path.append( os.path.curdir )

     
     
#
# make the following defines in the main module
#
exec("""
#
# needed for command tracing
#
import trace     

#
# preload these container types
#
import vec3
import cdsList
import cdsVector
import atomSel
#import xplorPot
import atom
from atomSel import AtomSel
from potList import PotList
from energyReport import EnergyReport

import simulationWorld
simWorld = simulationWorld.SimulationWorld_world()

#augment the help function a bit
import xplorDoc

help_xplorDoc_dontCall=1
def help(obj=0):
  if obj:
    xplorDoc.help(obj)
  else:
    xplorDoc.help()

""",__main__.__dict__)

#
# add convenience access to xplor.wrap.en/disableOutput()
#
from xplor import wrap
enableOutput=wrap.enableOutput
disableOutput=wrap.disableOutput

#
# exit message
#
exitMessage=''
def writeExitMessage():
    import simulationWorld
    simWorld = simulationWorld.SimulationWorld_world()
    if simWorld.logLevel()!='none' and exitMessage:
        writeConsole(exitMessage)
        pass
    return
import atexit
atexit.register(writeExitMessage)

#
#
#
def writeConsole(msg):
    """write a message to the user's console - file descriptor number 9.
    This is stdout from where the xplor command is typed.
    """
    import os
    try:
        os.write(9,bytes(msg,"utf8"))
    except OSError:
        # fd 9 not open - this is the case for an extension
        pass
    

#
# insert parallel processing info, start up socket communications
#
from os import environ as env
p_processID = -1
p_numProcs = -1
p_host0 = env['XPLOR_PHOST']
p_port  = int(env['XPLOR_PPORT'])
allowStartupFailures=0
if 'XPLOR_ALLOWSTARTUPFAILURES' in env:
   allowStartupFailures = int(env['XPLOR_ALLOWSTARTUPFAILURES'])
from socketComm import Comm
numProcs=int( env['XPLOR_NUM_PROCESSES'] )
p_comm = Comm(numProcs,
              int( env['XPLOR_PROCESS'] ),p_host0,p_port,
              startupTimeout=50+numProcs*2,
              allowStartupFailures=allowStartupFailures,
              startupDelay=float( env['XPLOR_STARTUPDELAY'] ) *
                           int( env['XPLOR_PROCESS'] )
              )
p_numProcs = p_comm.numProcs
p_processID = p_comm.procNum
env['XPLOR_NUM_PROCESSES'] = '%d' % p_numProcs
env['XPLOR_PROCESS'] = '%d' % p_processID





#
# add random number accessors to xplor.simulation for backwards compatibility
#

exec("""
import simulationWorld
def randomSeed():
  return simulationWorld.SimulationWorld_world().random.seed()
def setRandomSeed(n):
  simulationWorld.SimulationWorld_world().setRandomSeed(n)
  return
def uniformRandom():
  simulationWorld.SimulationWorld_world().random.uniform()
""",simulation.__dict__)

