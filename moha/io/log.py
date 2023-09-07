import atexit
from contextlib import contextmanager
import datetime
from functools import wraps
import getpass
import os
import io
import sys
import resource
import time
import urllib

import moha

__all__ = ['log', 'timer']


class Log(object):
    """Output file.

    Attributes
    ----------
    width : int
        Text width.

    Methods
    -------
    __init__(self, name, version, head_banner, foot_banner, timer)
        Initialize a ScreenLog object.

    __call__(self, *words)
        Print nicely formatted output to screen.

    output(self,target)
        Specify the output target.

    hline(self, char='~')
        Print a horizontal line.

    blank(self)
        Print a blank line.

    print_header(self)
        Print the header.

    print_footer(self)
        Print the footer.

    _print_basic_info(self)
        Print basic runtime info.
    """
    width = 80

    def __init__(self, name, version, head_banner, foot_banner, timer):
        """Initialize a ScreenLog object.

        Parameters
        ----------
        name : str
            Name of the package.

        version : str
            Version of the package.

        head_banner : str
            The first text to be printed on target.

        foot_banner : str
            The last text to be printed on target.

        timer : Timer
            A timer to keep track of CPU time taken by different parts of the
            package.
        """
        self.name = name
        self.version = version
        self.head_banner = head_banner
        self.foot_banner = foot_banner
        self.timer = timer
        self.target = io.StringIO()

        self._active = False


    def __call__(self, *words):
        """Print nicely formatted output to screen.

        All arguments are joined with a spaced and wrapped to a fixed line width before
        printing. The first ampersand is interpreted as the place of a hanging indent.
        """
        if not self._active:
            self.print_header()
        print('{0:s}'.format(*words), file=self.target)

    def output(self,target):
        """Specify the output target.

        Parameters
        ----------
        target : str
            The output target.
        """
        if isinstance(target,str):
            if target=='silent':
                self.target = io.StringIO()
            elif target=='screen':
                self.target = sys.stdout
            elif target.endswith('.log'):
                self.target = open(target,'w')
            else:
                raise ValueError('The format of output file must be .log')
        else:
            raise TypeError("Output target not supported")
        
    def hline(self, char='~'):
        """Print a horizontal line.

        Parameters
        ----------
        char : str
            The line consists a repetition of this character.
        """
        self(char*self.width)

    def blank(self):
        """Print a blank line."""
        print(file=self.target)

    def print_header(self):
        """Print the header."""
        if not self._active:
            self._active = True
            print(self.head_banner, file=self.target)
            self._print_basic_info()

    def print_footer(self):
        """Print the footer."""
        self.timer._stop('Total')
        self.timer.report(self)
        print(self.foot_banner, file=self.target)

    def _print_basic_info(self):
        """Print basic runtime info."""
        #Grab basic info
        rows = []
        rows.append('User:               &'+getpass.getuser())
        rows.append('Machine Info:       &'+' '.join(os.uname()))
        rows.append('Time:               &'+datetime.datetime.now().isoformat())
        rows.append('Python Version:     &'+sys.version.replace('\n', ''))
        rows.append('Current Dir:        &'+os.getcwd())
        rows.append('Command Line:       &'+''.join(sys.argv))

        self.blank()
        #Check for alignment
        for row in rows:
            pos = row.find(u'&')
            if pos == -1:
                lead = u''
                rest = row
            else:
                lead = row[:pos] + ' '
                rest = row[pos+1:]

            width = self.width - len(lead)

            #Break and print the line
            first = True
            while len(rest) > 0:
                if len(rest) > width:
                    pos = rest.rfind(' ', 0, width)
                    if pos == -1:
                        current = rest[:width]
                        rest = rest[width:]
                    else:
                        current = rest[:pos]
                        rest = rest[pos:].lstrip()
                else:
                    current = rest
                    rest = u''
                print(u'{0:s}{1:s}'.format(lead, current), file=self.target)
                if first:
                    lead = u' '*len(lead)
                    first = False
        self.blank()
        self.hline(char='=')
        self.blank()

class Timer(object):
    """Keep track of CPU time of a part of the program.

    This is the bare-bones timer. It does not make any distinction between time spent in
    its own code or code of other functions it calls (if they are also timed).
    """

    def __init__(self):
        """Initialize the Timer object."""
        self.cpu = 0.0
        self._start = None
        # The _depth attribute is needed for timed recursive functions.
        self._depth = 0

    def start(self):
        """Mark start of a function."""
        if self._depth == 0:
            assert self._start is None
            self._start = time.time()
        self._depth += 1

    def stop(self):
        """Mark end of a function."""
        if self._depth > 0:
            assert self._start is not None
            self._depth -= 1
        if self._depth == 0:
            self.cpu += time.time() - self._start
            self._start = None


class SubTimer(object):
    """Keep track of CPU time in a part of a program, aware of timed called functions."""

    def __init__(self, label):
        """Initialize a SubTimer object.

        Parameters
        ----------
        label : str
            A name for this timer.
        """
        self.label = label
        self.total = Timer()
        self.own = Timer()

    def start(self):
        """Mark start of the timed function."""
        self.total.start()
        self.own.start()

    def start_sub(self):
        """Mark start of a timed function inside the current function."""
        self.own.stop()

    def stop_sub(self):
        """Mark end of a timed function inside the current function."""
        self.own.start()

    def stop(self):
        """Mark end of the timed function."""
        self.own.stop()
        self.total.stop()


class TimerGroup(object):
    """Keep track of CPU time spent in different parts of a code."""

    def __init__(self):
        """Initialize a TimerGroup object."""
        self.reset()

    def reset(self):
        """Reset all timers."""
        self.parts = {}
        self._stack = []
        self._start('Total')

    @contextmanager
    def section(self, label):
        """Time part of a code with a contextmanager.

        Parameters
        ----------
        label : str
            A label for that part of the code.

        This can be used as follows:

        .. code-block:: python

            with log.timer('AlgoWhatever'):
                # Computationally expensive code here, instead of pass:
                pass
            # Back to code that is has no significant cost.
        """
        self._start(label)
        try:
            yield
        finally:
            self._stop(label)

    def with_section(self, label):
        """Time part of a code with a decorator.

        Parameters
        ----------
        label : str
            A label for that part of the code.

        This can be used as follows:

        .. code-block:: python

            @log.timer.with_section('AlgoThisAndThat')
            def slow_function(*args, **kwargs):
                # Slow code inside...
                pass
        """
        def decorator(fn):
            """Decorator that adds timer to function."""
            @wraps(fn)
            def wrapper(*args, **kwargs):
                """Wrap function with timer."""
                with self.section(label):
                    return fn(*args, **kwargs)
            return wrapper
        return decorator

    def _start(self, label):
        """Start timing part of a code.

        This method should not be called directly. Use decorator and contextmanager API
        instead.

        Parameters
        ----------
        label : str
            A label for that part of the code.
        """
        assert len(label) <= 14
        # get the right timer object
        timer = self.parts.get(label)
        if timer is None:
            timer = SubTimer(label)
            self.parts[label] = timer
        # start timing
        timer.start()
        if len(self._stack) > 0:
            self._stack[-1].start_sub()
        # put it on the stack
        self._stack.append(timer)

    def _stop(self, label):
        """Stop timing part of a code.

        This method should not be called directly. Use decorator and contextmanager API
        instead.

        Parameters
        ----------
        label : str
            A label for that part of the code.
        """
        timer = self._stack.pop(-1)
        assert timer.label == label
        timer.stop()
        if len(self._stack) > 0:
            self._stack[-1].stop_sub()

    def report(self, log):
        """Write a report of the CPU usage on screen.

        Parameters
        ----------
        log : ScreenLog
            A logger to use for writing the screen output.
        """
        log.hline(char='=')
        log.blank()
        log('Overview of CPU time usage.')
        log.hline()
        log('{0:<20s}{1:<20s}{2:<20s}'.format('Label','Total','Own'))
        log.hline()
        for label, timer in sorted(self.parts.items()):
            log('{0:20s} {1:<20.1f} {2:<20.1f}'.format(
                label,
                timer.total.cpu, 
                timer.own.cpu,
            ))
        log.hline()
        ru = resource.getrusage(resource.RUSAGE_SELF)
        log('{0:<20s}{1:<20.1f}'.format('CPU user time', ru.ru_utime))
        log('{0:<20s}{1:<20.1f}'.format('CPU system time', ru.ru_stime))
        log('{0:<20s}{1:<20.1f}'.format('Page swaps', ru.ru_nswap))



head_banner = """\
================================================================================

Welcome to MOHA %s!

MOHA is written and maintained by Yilin Zhao (1).

(1) The Ayers Group, McMaster University, Hamilton, Ontario, Canada.

More information about MOHA can be found on this website:
https://zhaoyilin.github.io/moha/

The purpose of this log file is to track the progress and quality of a 
computation. 

================================================================================""" % (moha.__version__)


foot_banner = """
================================================================================

End of the MOHA program.

Thank you for using MoHa %s!

================================================================================""" % (moha.__version__)

timer = TimerGroup()
log = Log('MOHA', 'moha.__version__', head_banner, foot_banner, timer)

atexit.register(log.print_footer)





