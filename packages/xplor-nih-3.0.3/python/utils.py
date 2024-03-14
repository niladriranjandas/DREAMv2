"""Module with miscellaneous utilities.


"""

import datetime
import calendar


# Key: month name (e.g., 'Jan', 'Apr'), value: the month numerical value.
moname2num = dict((v,k) for k,v in enumerate(calendar.month_abbr))


def duration(start, end):
    """Return duration (in seconds) between start and end times.

    This function is useful for timing the duration of an Xplor-NIH run.  The
    arguments start and end are strings specifying the start and end times of
    the run.  Their format is of the type 'hh:mm:ss dd-Mmm-yy', e.g., as in
    '14:07:58 17-Jun-15' (this is the format used by the "entry time" or "exit
    time" at end of an Xplor-NIH logfile).

    It is assumed the years are in the 2000's.

    """
    # Start.
    date = start.split()[1].split('-')
    time = start.split()[0].split(':')
    
    (year, month, day) = (2000+int(date[2]), moname2num[date[1]], int(date[0]))
    (hours, minutes, seconds) = (int(time[0]), int(time[1]), int(time[2]))

    start = datetime.datetime(year, month, day, hours, minutes, seconds)

    # End.
    date = end.split()[1].split('-')
    time = end.split()[0].split(':')
    
    (year, month, day) = (2000+int(date[2]), moname2num[date[1]], int(date[0]))
    (hours, minutes, seconds) = (int(time[0]), int(time[1]), int(time[2]))

    end = datetime.datetime(year, month, day, hours, minutes, seconds)

    duration = end - start
    return duration.total_seconds()

                          
def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2)+1):
        yield chr(c)
        pass
    pass

def printReStructuredText(rstString,
                          outputType="txt",
                          cachePrefix=None,
                          cacheDir=None,
                          cacheTime=None):
    """Print out reStructuredText as raw rst, markup-stripped txt or displayed
    as a pdf, depending on the outputType argument. The cache- arguments
    configured location and names for a cached PDF version which is regenerated
    if cacheTime > modification time of the cached PDF.
    """

    if outputType=="rst":
        print(rstString)
    elif outputType=="txt":
        from docutils.core import publish_string
        import rst2txt
        print(publish_string(rstString, writer=rst2txt.Writer(),
                             settings_overrides={'output_encoding': 'unicode'})
              )
    elif outputType=="pdf":
        
        import os
        if not cacheDir or not os.path.exists(cacheDir):
            print("printReStructuredText: Error. no such directory:", end=' ')
            print(cacheDir)
            pass
        
        os.chdir(cacheDir)

        cachePrefix=cachePrefix.strip()

        pdfFilename=cachePrefix + '.pdf'

        regen = True
        if os.path.exists(pdfFilename):
            pdfModTime = os.stat(pdfFilename).st_mtime
            if cacheTime < pdfModTime:
                regen = False
                pass
            pass
        if regen:
            try:
                from docutils.core import publish_string
                from docutils.writers.latex2e import Writer
                texFilename = cachePrefix + '.tex'
                open(texFilename,'w').write(
                    publish_string(rstString,
                                   settings_overrides={
                    'documentoptions':'12pt',
#                    'stylesheet':'stylesheet',
                    'output_encoding': 'unicode'
                    },
                                   writer=Writer()))
                import subprocess
                try:
                    cmd=["pdflatex",texFilename]
                    subprocess.run(cmd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   timeout=10)
                except subprocess.SubprocessError as err:
                    print("ERROR running: "+" ".join(cmd))
                    print("output:", err.stdout.decode("utf-8"))
                    os.unlink(pdfFilename)
                    pass
                pass
            except KeyError:
                print("Warning: could not regenerate PDF for: ", cachePrefix)
                pass
            pass
        try:
            from os import environ as env
            openCommand = env["PDF_OPEN_COMMAND"]
        except:
            print("Error reading environment variable PDF_OPEN_COMMAND")
            raise
        os.system(openCommand + " " +pdfFilename)
        #FIX: change back to original directory
    else:
        print("printReStructuredText: Error. bad outputType:", outputType)
        pass
    return

        
        
