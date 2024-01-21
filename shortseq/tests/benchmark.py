import os
import gzip
import inspect
import unittest
import numpy as np
import matplotlib as mpl
default_backend = mpl.get_backend()
import matplotlib.pyplot as plt

# umi_tools changes backend upon import
from umi_tools._dedup_umi import edit_distance
mpl.use(default_backend)

from scipy.spatial.distance import hamming
from typing import Dict, List, Iterable
from pympler.asizeof import asizeof
from collections import defaultdict
from datetime import datetime
from random import randint
from timeit import timeit
from glob import glob

import shortseq as sq
from util import rand_sequence, sorted_natural

import warnings
warnings.filterwarnings("error")
np.seterr(all='raise')

# === TYPE HINTS =============================================================

DataSet = Dict[int, List[bytes]]
FinalResults = Dict[str, Iterable[float]]


# === BENCHMARKS ==============================================================

class MemoryBenchmarks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = make_data()

    def test_mem_by_length(self):
        """Memory usage by sequence length is measured using asizeof() for
        each object type. NOTE:
            1) The measurement of Gzip Bytes is the number of bytes in
               the compressed sequence, NOT the size of the PyBytes
               object that gzip.compress() returns.
            2) The NumPy array is constructed as a tokenized form of the
               sequence that is amenable to edit distance calculation.
        """

        title = "Memory Usage by Object Type"
        lab_x = "Sequence Length"
        lab_y = "Object Size (bytes)"

        mem_sq = {}
        mem_np = {}
        mem_un = {}
        mem_by = {}
        mem_gz = {}

        for length, seqs in self.data.items():
            mem_sq[length] = asizeof(sq.pack(seqs[0]))
            mem_np[length] = asizeof(np.char.asarray(seqs[0], itemsize=1))
            mem_un[length] = asizeof(str(seqs[0]))
            mem_by[length] = asizeof(seqs[0])
            mem_gz[length] = len(gzip.compress(seqs[0]))

        results = {
            'ShortSeq':   mem_sq.values(),
            'NumPy':      mem_np.values(),
            'PyUnicode':  mem_un.values(),
            'PyBytes':    mem_by.values(),
            'Gzip Bytes': mem_gz.values(),
        }

        save_and_plot(results, title, lab_x, lab_y)


class TimeBenchmarks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = make_data()

    def test_construction_from_bytes(self):
        """Measures the amount of time that it takes to construct a sequence object from a PyBytes object.
        PyBytes is the natural choice over PyUnicode for the initial sequence object because it is more
        memory efficient and faster to construct when reading from a file (mode "rb"). In order to be
        consistent with the edit distance benchmark, NumPy's construction time is measured for a
        tokenized form that is amenable to edit distance calculation."""

        samples = 1000
        title = "Construction from PyBytes"
        lab_x = "Sequence Length"
        lab_y = "Average Time (10$^{-6}$ seconds)"

        times_sq = defaultdict(list)
        times_py = defaultdict(list)
        times_dc = defaultdict(list)
        times_np = defaultdict(list)

        for length, seqs in self.data.items():
            for seq in seqs:
                times_sq[length].append(timeit(lambda: sq.pack(seq), number=samples) / samples)
                times_py[length].append(timeit(lambda: str(seq), number=samples) / samples)
                times_dc[length].append(timeit(lambda: seq.decode(), number=samples) / samples)
                times_np[length].append(timeit(lambda: np.char.asarray(seq, itemsize=1), number=samples) / samples)

        results = {
            'sq.pack(x)':         [sum(t) / len(t) for t in times_sq.values()],
            'str(x)':             [sum(t) / len(t) for t in times_py.values()],
            'x.decode()':         [sum(t) / len(t) for t in times_dc.values()],
            'np.char.asarray(x)': [sum(t) / len(t) for t in times_np.values()],
        }

        fig, ax = save_and_plot(results, title, lab_x, lab_y)

        # Validate y-axis label
        y_axis_mag = ax.yaxis.get_major_formatter().orderOfMagnitude
        self.assertEqual(y_axis_mag, -6, f"Expected y-axis order of magnitude to be -6, got {y_axis_mag} instead.")

    def test_hamming_distance(self):
        samples = 100
        title = "Edit Distance Calculation"
        lab_x = "Sequence Length"
        lab_y = "Average Time (seconds)"

        times_sq = defaultdict(list)
        times_np = defaultdict(list)
        times_um = defaultdict(list)
        times_py = defaultdict(list)
        times_sy = defaultdict(list)

        for length, seqs in self.data.items():
            indices = [(randint(0, len(seqs) - 1), randint(0, len(seqs) - 1)) for _ in range(10)]
            sq_seqs = [sq.pack(seq) for seq in seqs]
            np_seqs = [np.char.asarray(seq, itemsize=1) for seq in seqs]
            py_seqs = [str(seq) for seq in seqs]
            sy_seqs = [list(seq) or [""] for seq in seqs]

            for i, j in indices:
                a, b = sq_seqs[i], sq_seqs[j]
                c, d = py_seqs[i], py_seqs[j]
                e, f = sy_seqs[i], sy_seqs[j]
                g, h = seqs[i], seqs[j]
                k, l = np_seqs[i], np_seqs[j]

                times_sq[length].append(timeit(lambda: a ^ b, number=samples) / samples)
                times_np[length].append(timeit(lambda: np.count_nonzero(k != l), number=samples) / samples)
                times_um[length].append(timeit(lambda: edit_distance(g, h), number=samples) / samples)
                times_py[length].append(timeit(lambda: sum(a_nt != b_nt for a_nt, b_nt in zip(c, d)), number=samples) / samples)
                times_sy[length].append(timeit(lambda: hamming(e, f) * len(e), number=samples) / samples)

        results = {
            'ShortSeq':  [sum(t) / len(t) for t in times_sq.values()],
            'NumPy':     [sum(t) / len(t) for t in times_np.values()],
            'UMI-tools': [sum(t) / len(t) for t in times_um.values()],
            'PyUnicode': [sum(t) / len(t) for t in times_py.values()],
            'SciPy':     [sum(t) / len(t) for t in times_sy.values()],
        }

        save_and_plot(results, title, lab_x, lab_y, log_scale=True)


class ReplotTests(unittest.TestCase):
    ...

    # def test_template_for_replotting(self):
    #     # Specify benchmark data
    #     test_name = ""
    #     timestamp = ""  # optional, defaults to most recent run
    #
    #     # Plot attributes
    #     title = ""
    #     lab_x = ""
    #     lab_y = ""
    #
    #     prev, run_dir = load_results(test_name, timestamp)
    #     fig, ax = plot_results(prev, title, lab_x, lab_y, outfile=None)
    #
    #     # Modify plot figure or axes
    #     ...
    #
    #     # Save plot
    #     ts = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    #     fig.savefig(os.path.join(run_dir, f"plot_{ts}.svg"))


# === HELPER FUNCTIONS ========================================================


def make_data(min_len=0, max_len=1024, n_samples=10) -> DataSet:
    """Generate a dictionary of random sequences of varying lengths, with n_samples per length."""

    return {
        length: [
            rand_sequence(length, as_bytes=True)
            for _ in range(n_samples)
        ]
        for length in range(min_len, max_len+1)
    }


def save_and_plot(sets: FinalResults, title, x_label, y_label, dirname=None, **kwargs):
    """Saves the benchmark results and plot to a timestamped subdirectory within benchmarks/{dirname}."""

    if dirname is None:  # just use the test function's name
        dirname = inspect.stack()[1].function.replace("test_", "")

    subdir = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    outpath = os.path.join("benchmarks", dirname, subdir)
    outfile = os.path.join(outpath, "plot.svg")
    os.makedirs(outpath, exist_ok=True)

    save_results(sets, outpath)
    return plot_results(sets, title, x_label, y_label, outfile=outfile, **kwargs)


def save_results(sets: FinalResults, outdir):
    """Saves the benchmark results as .txt files in outdir to allow for re-plotting."""

    for label, data in sets.items():
        with open(os.path.join(outdir, f"{label}.txt"), 'w') as f:
            f.write("\n".join(map(str, data)))


def plot_results(sets: FinalResults, title: str, x_label: str, y_label: str, **kwargs):
    """Plots a dictionary of data sets. Keys are legend labels, and values are lists of values to plot."""

    fig, ax = plt.subplots()
    ax: plt.Axes; fig: plt.Figure
    plt.rcParams['figure.dpi'] = 300
    if kwargs.get('log_scale', False):
        ax.tick_params(which="both", width=1.5, axis='y')
        ax.tick_params(which="major", length=4, axis='y')
        ax.set_yscale('log')

    for label, values in sets.items():
        ax.plot(values, label=label)

    plt.legend(loc='upper left')
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if kwargs.get('caption') is not None:
        cap = kwargs['caption']
        fig.text(0, -0.05, cap, ha='left', wrap=True)

    if kwargs.get('outfile') is not None:
        plt.savefig(kwargs['outfile'])

    plt.tight_layout()
    plt.show()  # must be called after savefig else a blank plot is saved
    return fig, ax


def load_results(test_name, timestamp=None):
    test_root = os.path.join("benchmarks", test_name)
    if not timestamp:
        timestamp = os.path.basename(sorted_natural(glob(f"{test_root}/*"))[-1])  # most recent run

    run_dir = os.path.join(test_root, timestamp)
    results = {}

    for fname in glob(f"{run_dir}/*.txt"):
        with open(fname, 'r') as f:
            label = os.path.basename(fname).replace(".txt", "")
            values = [float(v) for v in f.readlines()]
            results[label] = values

    return results, run_dir


if __name__ == '__main__':
    unittest.main()
