import contextlib
import io
import itertools
import os
import pathlib
import pty
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam

from align_trim.console import Console, bam_fraction, console
from align_trim.main import go
from align_trim.qc import QC
from tests.test_integration import (
    BAM_PATH_PAIRED_V3_0_0,
    BAM_PATH_V5_3_2,
    BED_PATH_V3_0_0,
    create_args,
)


class FakeTTY(io.StringIO):
    encoding = "utf-8"

    def isatty(self):
        return True


def fake_clock(step):
    ticks = itertools.count()
    return lambda: next(ticks) * step


class ConsoleTestCase(unittest.TestCase):
    def setUp(self):
        env = mock.patch.dict(os.environ)
        env.start()
        self.addCleanup(env.stop)
        for key in ("NO_COLOR", "FORCE_COLOR"):
            os.environ.pop(key, None)
        os.environ["TERM"] = "xterm-256color"

    def capture(self, fn, stream=None):
        stream = stream if stream is not None else io.StringIO()
        with contextlib.redirect_stderr(stream):
            fn()
        return stream.getvalue()


class TestConsole(ConsoleTestCase):
    def test_plain_output_in_a_log(self):
        c = Console()

        def write():
            c.heading("Coverage")
            c.prefixed("QC:", "1 fragments")
            c.warn("x")
            c.error("bad, exiting.")

        out = self.capture(write)
        self.assertEqual(
            out, "== Coverage ==\nQC: 1 fragments\nWARNING: x\nbad, exiting.\n"
        )
        self.assertNotIn("\x1b", out)
        self.assertNotIn("\r", out)

    def test_colour_on_a_terminal(self):
        c = Console()
        out = self.capture(lambda: (c.heading("Coverage"), c.warn("x")), FakeTTY())
        self.assertIn("\x1b[", out)
        self.assertIn("─", out)
        self.assertIn("WARNING:", out)

    def test_no_color(self):
        os.environ["NO_COLOR"] = "1"
        c = Console()
        out = self.capture(lambda: c.warn("x"), FakeTTY())
        self.assertEqual(out, "WARNING: x\n")

    def test_force_color(self):
        os.environ["FORCE_COLOR"] = "1"
        c = Console()
        out = self.capture(lambda: (c.heading("Run"), c.warn("x")))
        self.assertIn("\x1b[", out)
        # colour is forced but a log still gets plain headings
        self.assertIn("== Run ==", out)

    def test_dumb_terminal(self):
        os.environ["TERM"] = "dumb"
        c = Console()
        out = self.capture(lambda: (c.heading("Run"), c.warn("x")), FakeTTY())
        self.assertEqual(out, "== Run ==\nWARNING: x\n")

    def test_heading_not_repeated(self):
        c = Console()
        out = self.capture(lambda: (c.heading("Run"), c.heading("Run")))
        self.assertEqual(out, "== Run ==\n")

    def test_warnings_flushed_under_heading_once(self):
        c = Console()
        self.assertEqual(self.capture(c.flush_warnings), "")
        c.defer_warning("one")
        c.read_warning("read1: problem soft masking")
        out = self.capture(c.flush_warnings)
        self.assertEqual(
            out,
            "== Warnings ==\n"
            "WARNING: 1 fragments skipped because primer soft masking failed; rerun with --verbose for details\n"
            "WARNING: one\n",
        )
        self.assertEqual(self.capture(c.flush_warnings), "")

    def test_read_warnings_only_printed_when_verbose(self):
        c = Console()
        self.assertEqual(self.capture(lambda: c.read_warning("read1: x")), "")
        c.reset(verbose=True)
        self.assertEqual(
            self.capture(lambda: c.read_warning("read1: x")), "warning  read1: x\n"
        )


class TestProgress(ConsoleTestCase):
    def test_live_bar_on_a_terminal(self):
        c = Console()
        stream = FakeTTY()

        def run():
            p = c.progress(fraction=lambda: 0.5, clock=fake_clock(0.05))
            self.assertEqual(sum(1 for _ in p.track(range(100_000))), 100_000)
            p.done()

        out = self.capture(run, stream)
        redraws = out.count("trimming")
        # checked every 8192 records, redrawn at most every 0.125s
        self.assertGreaterEqual(redraws, 1)
        self.assertLessEqual(redraws, 100_000 // 8192)
        self.assertIn("50%", out)
        self.assertIn("ETA", out)
        self.assertEqual(out.count("processed 100,000 records"), 1)
        # the bar is cleared before the completion line
        before = out.split("processed")[0]
        self.assertIn("\r\x1b[K", before[-20:])
        self.assertIn("✔", before[-20:])

    def test_count_without_a_total(self):
        c = Console()
        out = self.capture(
            lambda: list(
                c.progress(fraction=None, clock=fake_clock(1)).track(range(20_000))
            ),
            FakeTTY(),
        )
        self.assertIn("records", out)
        self.assertIn("elapsed", out)
        self.assertNotIn("ETA", out)

    def test_periodic_lines_in_a_log(self):
        c = Console()

        def run():
            p = c.progress(fraction=lambda: 0.25, clock=fake_clock(10))
            list(p.track(range(100_000)))
            p.done("5 read pairs")

        out = self.capture(run)
        self.assertNotIn("\r", out)
        self.assertNotIn("\x1b", out)
        lines = out.splitlines()
        self.assertTrue(lines[0].startswith("progress: "), lines)
        self.assertIn("(25%)", lines[0])
        self.assertRegex(
            lines[-1], r"^processed 100,000 records \(5 read pairs\) in [\d.]+s, "
        )
        self.assertEqual(sum(1 for line in lines if line.startswith("processed")), 1)

    def test_verbose_disables_the_live_bar(self):
        c = Console()
        c.reset(verbose=True)
        out = self.capture(
            lambda: list(c.progress(clock=fake_clock(0.05)).track(range(20_000))),
            FakeTTY(),
        )
        self.assertNotIn("\r", out)

    def test_no_live_bar_when_sam_is_on_the_terminal(self):
        c = Console()
        out = self.capture(
            lambda: list(
                c.progress(live=False, clock=fake_clock(10)).track(range(100_000))
            ),
            FakeTTY(),
        )
        self.assertNotIn("\r", out)
        self.assertIn("progress: ", out)

    def test_quiet(self):
        c = Console()
        c.reset(quiet=True)

        def run():
            p = c.progress(clock=fake_clock(10))
            list(p.track(range(100_000)))
            p.done()
            self.assertEqual(p.count, 100_000)

        self.assertEqual(self.capture(run), "")

    def test_line_clears_the_live_bar(self):
        c = Console()

        def run():
            p = c.progress(clock=fake_clock(1))
            for i in p.track(range(20_000)):
                if i == 10_000:
                    c.line("message")

        out = self.capture(run, FakeTTY())
        self.assertIn("\r\x1b[Kmessage\n", out)

    def test_bam_fraction(self):
        with pysam.AlignmentFile(str(BAM_PATH_V5_3_2)) as infile:
            fraction = bam_fraction(infile, str(BAM_PATH_V5_3_2))
            self.assertLess(fraction(), 0.5)
            for _ in infile:
                pass
            self.assertGreaterEqual(fraction(), 0.95)
            self.assertIsNone(bam_fraction(infile, "-"))
            self.assertIsNone(bam_fraction(infile, None))


class TestSections(ConsoleTestCase):
    def run_go(self, **kwargs):
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory(dir="tests", suffix="-console") as tempdir:
            args = create_args(
                output=(pathlib.Path(tempdir) / "out.sam").absolute(),
                qc_prefix=f"{tempdir}/qc",
                genome_coverage_report=f"{tempdir}/cov",
                bedfile=BED_PATH_V3_0_0,
                samfile=BAM_PATH_PAIRED_V3_0_0,
                normalise=200,
                **kwargs,
            )
            write_overlaps = QC._write_overlaps

            # the test data has no contamination, so add a QC warning
            def with_warning(qc, path):
                write_overlaps(qc, path)
                console.defer_warning(
                    "majority base differs between overlapping amplicons test"
                )

            with mock.patch.object(QC, "_write_overlaps", with_warning):
                with contextlib.redirect_stderr(stderr):
                    go(args)
        return stderr.getvalue()

    def headings(self, out):
        return [line for line in out.splitlines() if line.startswith("== ")]

    def test_section_order(self):
        out = self.run_go()
        self.assertEqual(
            self.headings(out),
            ["== Run ==", "== Coverage ==", "== QC summary ==", "== Warnings =="],
        )
        self.assertNotIn("\x1b", out)
        self.assertNotIn("\r", out)
        self.assertIn("proper pairs) in ", out)
        # existing lines are unchanged under the headings
        self.assertIn("Pre-normalisation:   >= 1x: ", out)
        self.assertIn("Post-normalisation: Genome coverage summary", out)
        self.assertIn("QC: split reads: ", out)
        self.assertIn(
            "WARNING: majority base differs between overlapping amplicons test", out
        )

    def test_quiet_keeps_summaries(self):
        out = self.run_go(quiet=True)
        self.assertEqual(
            self.headings(out),
            ["== Coverage ==", "== QC summary ==", "== Warnings =="],
        )
        self.assertNotIn("processed", out)


@unittest.skipUnless(hasattr(os, "openpty"), "needs a pty")
class TestStdoutUntouched(ConsoleTestCase):
    def test_sam_on_stdout_has_no_escape_codes(self):
        master, slave = pty.openpty()
        with tempfile.TemporaryFile() as stdout:
            proc = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    "from align_trim.main import main; main()",
                    str(BED_PATH_V3_0_0),
                    "--samfile",
                    str(BAM_PATH_PAIRED_V3_0_0),
                    "--normalise",
                    "50",
                ],
                stdout=stdout,
                stderr=slave,
                env=dict(os.environ, FORCE_COLOR="1"),
            )
            os.close(slave)
            # drain the pty while the process runs so it never blocks on a full buffer
            stderr = b""
            while True:
                try:
                    chunk = os.read(master, 65536)
                except OSError:
                    break
                if not chunk:
                    break
                stderr += chunk
            os.close(master)
            self.assertEqual(proc.wait(timeout=60), 0, stderr)
            stdout.seek(0)
            sam = stdout.read()

        self.assertNotIn(b"\x1b", sam)
        self.assertNotIn(b"\r", sam)
        self.assertTrue(sam.startswith(b"@"))
        self.assertIn(b"\x1b[", stderr)
        self.assertIn("── Run".encode(), stderr)


if __name__ == "__main__":
    unittest.main()
