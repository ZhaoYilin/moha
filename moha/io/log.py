# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Screen logging, timing and citation management

The goal of the screen logger is to track the progress of a computation in a convenient
human-readable way, possibly highlighting problematic situations. It is not intended as a
computer-readable output file that contains all the results of a computation. For that
purpose, all useful information is written to a binary checkpoint file or kept in memory
as attributes of the HORTON objects.
"""

import atexit
from contextlib import contextmanager
import datetime
from functools import wraps
import getpass
import os
import sys
import resource
import time
import urllib

import moha

__all__ = ['log', 'timer']


class ScreenLog(object):
    """Nice and configurable screen output.

    The screen logger can be configured to print more or less output, with a parameter
    called the *log level*. This can be one of the class attributes:

    - `silent`: nothing at all
    - `warning`: only report things that may be going wrong
    - `low`: terse screen output
    - `medium`: normal screen output (default)
    - `high`: plenty of detail
    - `debug`: too much detail for a regular user

    Output is printed by `ScreenLog` with its `__call__` method (and also a few others).
    The program that uses `ScreenLog` is responsible for checking the log level before
    generating screen output with the `do_*` properties.
    """

    # Log levels.
    silent = 0
    warning = 1
    low = 2
    medium = 3
    high = 4
    debug = 5

    # Screen width parameter.
    width = 100

    def __init__(self, name, version, head_banner, foot_banner, timer, f=None):
        """Initialize a ScreenLog object.

        Parameters
        ----------
        name : str
            Name of the program.
        version : str
            Version number of the program.
        head_banner : str
            The first text to be printed on screen.
        foot_banner : str
            The last test to be printed on screen.
        timer : Timer
            A timer to keep track of CPU time taken by different parts of the program.
        f : file
            When given, print log to this file instead of stdout.
        """
        self.name = name
        self.version = version
        self.head_banner = head_banner
        self.foot_banner = foot_banner
        self.timer = timer

        self._active = False
        self._level = self.medium
        self._last_blank = False
        self.add_newline = False
        if f is None:
            _file = sys.stdout
        else:
            _file = f
        self._file = _file

    # The following properties can be used by the program to see if output of a given
    # priority should be printed.
    do_warning = property(lambda self: self._level >= self.warning)
    do_low = property(lambda self: self._level >= self.low)
    do_medium = property(lambda self: self._level >= self.medium)
    do_high = property(lambda self: self._level >= self.high)
    do_debug = property(lambda self: self._level >= self.debug)

    def set_level(self, level):
        """Set the log level.

        Parameters
        ----------
        level : int
            The desired level of output.
        """
        if level < self.silent or level > self.debug:
            raise ValueError('The level must be one of the ScreenLog attributes.')
        self._level = level

    def with_level(self, level):
        """Decorator to change the log level inside a specific function.

        Parameters
        ----------
        level : int
            The desired level of output.
        """
        def decorator(fn):
            """Return a wrapped function in which the log level is set to the desired value."""
            @wraps(fn)
            def wrapper(*args, **kwargs):
                """The wrapper."""
                old_level = self._level
                self.set_level(level)
                try:
                    result = fn(*args, **kwargs)
                finally:
                    self.set_level(old_level)
                return result
            return wrapper
        return decorator

    def __call__(self, *words):
        """Print nicely formatted output to screen.

        All arguments are joined with a spaced and wrapped to a fixed line width before
        printing. The first ampersand is interpreted as the place of a hanging indent.
        """
        s = ' '.join(w for w in words)
        if not self.do_warning:
            raise RuntimeError('The runlevel should be at least warning when logging.')
        if not self._active:
            self.print_header()

        # Check for alignment code '&'
        pos = s.find(u'&')
        if pos == -1:
            lead = u''
            rest = s
        else:
            lead = s[:pos] + ' '
            rest = s[pos+1:]
        width = self.width - len(lead)
        if width < self.width/2:
            raise ValueError('The lead may not exceed half the width of the terminal.')

        # Break and print the line
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
            print(u'%s%s' % (lead, current), file=self._file)
            if first:
                lead = u' '*len(lead)
                first = False

        self._last_blank = False

    def warn(self, *words):
        """Format and print a warning.

        See `__call__` for more details.
        """
        self.blank()
        text = '!WARNING!&'+' '.join(w for w in words)
        self(text)
        self.blank()

    def hline(self, char='~'):
        """Print a horizontal line.

        Parameters
        ----------
        char : str
            The line consists a repetition of this character.
        """
        self(char*self.width)
    
    def center(self, *words, **kwargs):
        """Print centered text, optionally with an edge surrounding the text.

        See __call__ for the details. One keyword argument `edge` is accepted that will be
        printed to the left and the right of the centered text.
        """
        if len(kwargs) == 0:
            edge = ''
        elif len(kwargs) == 1:
            if 'edge' not in kwargs:
                raise TypeError('Only one keyword argument is allowed, that is edge')
            edge = kwargs['edge']
        else:
            raise TypeError('Too many keyword arguments. Should be at most one.')
        s = ' '.join(w for w in words)
        if len(s) + 2*len(edge) > self.width:
            raise ValueError('Line too long. center method does not support wrapping.')
        self('%s%s%s' % (edge, s.center(self.width-2*len(edge)), edge))

    def blank(self):
        """Print a blank line."""
        if not self._last_blank:

            print(file=self._file)
            self._last_blank = True

    def deflist(self, l):
        """Print a definition list.

        Parameters
        ----------
        l : list
            A list of keyword and value pairs. A table will be printed where the first
            column contains the keywords and the second column contains (wrapped) values.
        """
        widest = max(len(item[0]) for item in l)
        for name, value in l:
            self('  %s :&%s' % (name.ljust(widest), value))

    def progress(self, niter):
        """Return a progress bar for `niter` iterations.

        The progress bar can be used as follows:

        .. code-block:: python

            pb = log.progress(100)
            for i in xrange(100):
                pb()

        A progress bar is only active at the medium level.
        """
        return ProgressBar(niter, self._file, self.width, silent=self._level != self.medium)

    def print_header(self):
        """Print the first screen output."""
        # Suppress any logging as soon as an exception is not caught.
        def excepthook_wrapper(type, value, traceback):
            """Silence the logger (as soon as an exception is raised)."""
            self.set_level(self.silent)
            sys.__excepthook__(type, value, traceback)
        sys.excepthook = excepthook_wrapper

        if self.do_warning and not self._active:
            self._active = True

            print(self.head_banner, file=self._file)
            self._print_basic_info()

    def print_footer(self):
        """Print the last screen output."""
        if self.do_warning and self._active:
            self._print_basic_info()
            self.timer._stop('Total')
            self.timer.report(self)
            print(self.foot_banner, file=self._file)

    def _print_basic_info(self):
        """Print basic runtime info."""
        if self.do_low:
            self.blank()
            self('User:          &' + getpass.getuser())
            self('Machine info:  &' + ' '.join(os.uname()))
            self('Time:          &' + datetime.datetime.now().isoformat())
            self('Python version:&' + sys.version.replace('\n', ''))
            self('Current Dir:   &' + os.getcwd())
            self('Command line:  &' + ' '.join(sys.argv))
            self.blank()
            self.hline(char='=')


class ProgressBar(object):
    """Simple progress bar for the screen logger."""

    def __init__(self, niter, f, width, silent):
        """Initialize a ProgressBar object.

        Parameters
        ----------
        niter : int
            The number of iterations to reach completion of something iterative.
        f : file object
            The file to write the progress bar to.
        width : int
            The screen width.
        silent : bool
            When `True`, nothing gets printed.
        """
        self.niter = niter
        self.f = f
        self.width = width
        self.silent = silent
        self.count = 0
        self.nchar = 0

    def __call__(self, inc=1):
        """Increment the progress bar counter."""
        self.count += inc
        if not self.silent:
            new_nchar = (self.count*self.width)/self.niter
            if new_nchar > self.nchar:
                self.f.write('>'*(new_nchar - self.nchar))
                self.f.flush()
                self.nchar = new_nchar
            if self.count == self.niter:
                self.f.write('\n')
        elif self.count > self.niter:
            raise ValueError('Progress bar overflow.')


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

    def get_max_own_cpu(self):
        """Return the CPU time of the slowest part of the code."""
        result = None
        for part in self.parts.values():
            if result is None or result < part.own.cpu:
                result = part.own.cpu
        return result

    def report(self, log):
        """Write a report of the CPU usage on screen.

        Parameters
        ----------
        log : ScreenLog
            A logger to use for writing the screen output.
        """
        max_own_cpu = self.get_max_own_cpu()
        #if max_own_cpu == 0.0:
        #    return
        log.blank()
        log('Overview of CPU time usage.')
        log.hline()
        log('Label             Total      Own')
        log.hline()
        bar_width = log.width-33
        for label, timer in sorted(self.parts.items()):
            if max_own_cpu > 0:
                cpu_bar = "W"*int(timer.own.cpu/max_own_cpu*bar_width)
            else:
                cpu_bar = ""
            log('%14s %8.1f %8.1f %s' % (
                label.ljust(14),
                timer.total.cpu, timer.own.cpu, cpu_bar.ljust(bar_width),
            ))
        log.hline()
        ru = resource.getrusage(resource.RUSAGE_SELF)
        log.deflist([
            ('CPU user time', '% 10.2f' % ru.ru_utime),
            ('CPU system time', '% 10.2f' % ru.ru_stime),
            ('Page swaps', '% 10i' % ru.ru_nswap),
        ])
        log.hline()



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
log = ScreenLog('MOHA', 'moha.__version__', head_banner, foot_banner, timer)



atexit.register(log.print_footer)





