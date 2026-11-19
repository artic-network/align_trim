"""Terminal output for align_trim: section headings, colour and progress on stderr.

Everything here writes to the current sys.stderr (looked up at write time, so
redirect_stderr in tests gets plain output). stdout is never touched. Colour and
the live progress line are only used when stderr is a terminal; log files get
plain text with `== Section ==` headings and no escape codes.
"""

import io
import os
import shutil
import sys
import time

BOLD = "1"
DIM = "2"
RED = "31"
GREEN = "32"
CYAN = "36"

# check the clock every 8192 records
CHECK_MASK = 0x1FFF
# minimum seconds between live redraws, and between progress lines in a log
REDRAW_INTERVAL = 0.125
LOG_INTERVAL = 30.0

SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
ASCII_SPINNER = "|/-\\"


def _human(n):
    for unit, size in (("M", 1e6), ("k", 1e3)):
        if n >= size:
            return f"{n / size:.2f}{unit}" if n < 10 * size else f"{n / size:.1f}{unit}"
    return str(int(n))


def _duration(seconds):
    seconds = int(seconds)
    hours, rest = divmod(seconds, 3600)
    minutes, seconds = divmod(rest, 60)
    return (
        f"{hours}:{minutes:02d}:{seconds:02d}" if hours else f"{minutes}:{seconds:02d}"
    )


def bam_fraction(infile, path):
    """Return a function giving how far through a BAM file on disk the reader is, or None.

    Uses the compressed byte offset (the high bits of the BGZF virtual offset) over
    the file size, so it works without an index. Not available for stdin or SAM.
    """
    if not path or path == "-" or not str(path).endswith(".bam"):
        return None
    try:
        size = os.path.getsize(path)
        infile.tell()
    except (OSError, ValueError, NotImplementedError):
        return None
    if not size:
        return None

    def fraction():
        try:
            return min(max((infile.tell() >> 16) / size, 0.0), 1.0)
        except (OSError, ValueError):
            return None

    return fraction


class Console:
    def __init__(self):
        self.reset()

    def reset(self, quiet=False, verbose=False):
        """Start a new run: clear warnings and counters."""
        self.quiet = quiet
        self.verbose = verbose
        self.warnings = []
        self.trim_failures = 0
        self._section = None
        self._live = None  # the Progress currently drawn on the terminal

    # terminal detection, evaluated against the current stderr
    @property
    def stream(self):
        return sys.stderr

    @property
    def is_tty(self):
        try:
            return self.stream.isatty()
        except (AttributeError, ValueError):
            return False

    @property
    def color(self):
        if os.environ.get("FORCE_COLOR"):
            return True
        return (
            self.is_tty
            and not os.environ.get("NO_COLOR")
            and os.environ.get("TERM") != "dumb"
        )

    @property
    def live(self):
        return self.is_tty and os.environ.get("TERM") != "dumb"

    @property
    def unicode(self):
        encoding = getattr(self.stream, "encoding", None) or ""
        return encoding.lower().replace("-", "") == "utf8"

    @property
    def width(self):
        """Columns of the terminal stderr is on (stdout is usually redirected), default 80."""
        try:
            columns = os.get_terminal_size(self.stream.fileno()).columns
        except (AttributeError, OSError, ValueError, io.UnsupportedOperation):
            columns = shutil.get_terminal_size((80, 24)).columns
        return columns or 80

    def style(self, text, *codes):
        if not self.color:
            return text
        return f"\x1b[{';'.join(codes)}m{text}\x1b[0m"

    # writing
    def _write(self, text):
        try:
            self.stream.write(text)
            self.stream.flush()
        except BrokenPipeError:
            pass

    def line(self, msg=""):
        """Write a line, clearing the live progress line first if one is drawn."""
        if self._live is not None:
            self._live.clear()
        self._write(f"{msg}\n")

    def heading(self, title):
        if title == self._section:
            return
        self._section = title
        if self.live:
            rule = "─" if self.unicode else "-"
            width = min(self.width, 80)
            text = f"{rule * 2} {title} "
            self.line(self.style(text + rule * max(width - len(text), 2), BOLD, CYAN))
        else:
            self.line(f"== {title} ==")

    def prefixed(self, prefix, msg):
        """A line with a dimmed prefix such as `QC:`, whose text is unchanged in plain mode."""
        self.line(f"{self.style(prefix, DIM)} {msg}")

    def warn(self, msg):
        self.line(f"{self.style('WARNING:', BOLD, RED)} {msg}")

    def error(self, msg):
        self.line(self.style(msg, BOLD, RED))

    def debug(self, msg):
        self.line(self.style(msg, DIM))

    def read_warning(self, msg):
        """A per-read problem: counted for a summary warning, printed only with --verbose."""
        self.trim_failures += 1
        if self.verbose:
            self.line(f"{self.style('warning', RED)}  {msg}")

    def defer_warning(self, msg):
        self.warnings.append(msg)

    def flush_warnings(self):
        if self.trim_failures:
            self.warnings.insert(
                0,
                f"{self.trim_failures} fragments skipped because primer soft masking failed; "
                "rerun with --verbose for details",
            )
            self.trim_failures = 0
        if self.warnings:
            self.heading("Warnings")
            for msg in self.warnings:
                self.warn(msg)
        self.warnings = []

    def progress(self, fraction=None, enabled=True, live=True, clock=time.monotonic):
        return Progress(
            self, fraction=fraction, enabled=enabled, live=live, clock=clock
        )


class Progress:
    """Count records from an iterator, drawing a live bar on a terminal.

    In a log (or with --verbose) it prints a plain progress line every LOG_INTERVAL
    seconds instead. `done()` prints a single completion line.
    """

    def __init__(
        self, console, fraction=None, enabled=True, live=True, clock=time.monotonic
    ):
        self.console = console
        self.fraction = fraction
        self.enabled = enabled and not console.quiet
        self.live = live and self.enabled and console.live and not console.verbose
        self.clock = clock
        self.count = 0
        self.frame = 0
        self.drawn = False
        self.start = clock()

    def track(self, iterable):
        if not self.enabled:
            for self.count, item in enumerate(iterable, 1):
                yield item
            return

        interval = REDRAW_INTERVAL if self.live else LOG_INTERVAL
        next_t = self.start + interval
        n = 0
        try:
            for item in iterable:
                n += 1
                if not n & CHECK_MASK:
                    t = self.clock()
                    if t >= next_t:
                        self.count = n
                        self._update(t)
                        next_t = t + interval
                yield item
        finally:
            self.count = n
            self.clear()

    def _rate(self, t):
        elapsed = t - self.start
        return self.count / elapsed if elapsed > 0 else 0.0

    def _update(self, t):
        if not self.live:
            fraction = self.fraction() if self.fraction else None
            pct = f" ({fraction:.0%})" if fraction is not None else ""
            self.console.line(
                f"progress: {self.count:,} records{pct}, {_human(self._rate(t))} records/s"
            )
            return
        self.console._live = self
        self._write(self._render(t))

    def _render(self, t):
        c = self.console
        frames = SPINNER if c.unicode else ASCII_SPINNER
        self.frame = (self.frame + 1) % len(frames)
        rate = self._rate(t)
        stats = f"{_human(self.count)} records  {_human(rate)} rec/s"
        fraction = self.fraction() if self.fraction else None
        if fraction is None:
            return f"{c.style(frames[self.frame], CYAN)} trimming  {stats}  {_duration(t - self.start)} elapsed"

        remaining = (t - self.start) * (1 - fraction) / fraction if fraction else 0
        tail = f" {fraction:4.0%}  {stats}  ETA {_duration(remaining)}"
        width = c.width
        bar_width = max(10, min(40, width - len(tail) - 14))
        filled = int(bar_width * fraction)
        full, empty = ("━", "─") if c.unicode else ("#", "-")
        bar = c.style(full * filled, GREEN) + c.style(empty * (bar_width - filled), DIM)
        return f"{c.style(frames[self.frame], CYAN)} trimming {bar}{tail}"

    def _write(self, text):
        self.console._write(f"\r\x1b[K{text}")
        self.drawn = True

    def clear(self):
        if self.drawn:
            self.console._write("\r\x1b[K")
            self.drawn = False
        if self.console._live is self:
            self.console._live = None

    def done(self, detail=None):
        """Print the completion line (not in quiet mode)."""
        self.clear()
        if self.console.quiet:
            return
        t = self.clock()
        c = self.console
        detail = f" ({detail})" if detail else ""
        tick = c.style("✔" if c.unicode else "*", BOLD, GREEN) + " " if c.live else ""
        c.line(
            f"{tick}processed {self.count:,} records{detail} in {t - self.start:.1f}s, "
            f"{_human(self._rate(t))} records/s"
        )


console = Console()
