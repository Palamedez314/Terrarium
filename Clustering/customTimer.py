############################################################################################
# Eigene Timer-Klasse 
############################################################################################

import time

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""
class CustomTimer:
    def __init__(self):
        self._start_times : dict[str,float] = {}
        self._timing_descriptions : dict[str,str] = {}
        self._cumul_times : dict[str,float] = {}

    def start_variable(self, varname:str, description:str=""):
        """Start timer with name varname"""
        assert type(varname) == str
        if varname in self._start_times.keys():
            raise TimerError(f"Timer is running. Use .stop()/.stop_variable to stop it")
        self._start_times[varname] = time.perf_counter()
        assert type(description) == str
        if description == "":
            description = varname
        self._timing_descriptions[varname] = description

    def print_variable(self, varname, val, cumul:bool=False):
        assert type(varname) == str
        if cumul:
            print(f"{self._timing_descriptions[varname]}: {val:0.4f} seconds (cumul)")
        else:
            print(f"{self._timing_descriptions[varname]}: {val:0.4f} seconds")

    def stop_variable(self, varname:str, print_single:bool=True, print_cumul:bool=False):
        """Stop timer with name varname"""
        assert type(varname) == str
        if varname not in self._start_times.keys():
            raise TimerError(f"Timer is not running. Use .start()/.start_variable() to start it")
        elapsed_time = time.perf_counter() - self._start_times.pop(varname)
        if print_single:
            self.print_variable(varname, elapsed_time, cumul=False)
        if print_cumul:
            cumul_time = self._cumul_times.pop(varname, 0) + elapsed_time
            self.print_variable(varname, cumul_time, cumul=False)
        else:
            self._cumul_times.pop(varname, 0)

    def print_cumul_variable(self, varname:str):
        if varname not in self._cumul_times.keys():
            raise TimerError(f"Timer has no recorded time yet. Use .start()/.start_variable() and .pause()/.pause_variable() to record time values")
        self.print_variable(varname, self._cumul_times[varname], cumul=True)

    def pause_variable(self, varname:str, print_single:bool=False, print_cumul:bool=False) :
        assert type(varname) == str
        if varname not in self._start_times.keys():
            raise TimerError(f"Timer is not running. Use .start()/.start_variable() to start it")
        elapsed_time = time.perf_counter() - self._start_times.pop(varname)
        if varname not in self._cumul_times.keys():
            self._cumul_times[varname] = elapsed_time
        else:
            self._cumul_times[varname] += elapsed_time
        if print_single:
            self.print_variable(varname, elapsed_time, cumul=False)
        if print_cumul:
            self.print_cumul_variable(varname)

    def start(self, description="Elapsed time"):
        """Start a new timer"""
        self.start_variable("standard", description)

    def stop(self, print_single:bool=True, print_cumul:bool=False):
        """Stop the timer, and report the elapsed time"""
        self.stop_variable("standard", print_single=print_single, print_cumul=print_cumul)

    def pause(self, print_single:bool=False, print_cumul:bool=False):
        self.pause_variable("standard", print_single=print_single, print_cumul=print_cumul)
    
    def print_cumul(self):
        self.print_cumul_variable("standard")