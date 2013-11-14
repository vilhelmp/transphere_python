

import os as _os
#~ import sys as _sys
#~ import subprocess as _subprocess
#~ import scipy as _scipy
#~ from matplotlib import pyplot as _plt





########################################################################
# GENERAL HELP FUNCTIONS (move to adavis_core)
def check_input(input_dictionary, input_defaults):
    return 0

def make_dirs(path):
    import errno
    # the makedirs function will raise a EEXIST error 
    # if the directory already exists, if so it returns False
    # if some other error is raise, it will raise that 
    try:
        _os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True


class ChangeDirectory:
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = _os.getcwd()
        _os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        _os.chdir(self.savedPath)

# Now you can enter the directory like this:

#~ with cd("~/Library"):
    #~ # we are in ~/Library
    #~ run some code
    #~ import subprocess
    #~ subprocess.call("ls")
