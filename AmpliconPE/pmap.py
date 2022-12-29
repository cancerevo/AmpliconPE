"""Parallel (multi-threaded) map function for python. 

Uses multiprocessing.Pool with error-resistant importing. There are two map
functions:

1) pmap(function, iterable) -> rapid fork-based multi-threaded map function.

2) low_memory_pmap(function, iterable) -> a more memory-efficient version
    intended for function calls that are individually long & memory-intensive.

"""
import os
import progressbar as pb
from time import sleep
import numpy as np
import pandas as pd
import traceback
import multiprocessing
from pathlib import Path
from multiprocessing.pool import Pool

CPUs = multiprocessing.cpu_count()


def simple_pmap(func, Iter, processes=CPUs - 1):
    with Pool(processes=processes) as P:
        return P.map(func, Iter)


def low_memory_pmap(func, Iter, processes=int(round(CPUs / 2)), chunksize=1):
    with Pool(processes=processes) as P:
        return [result for result in P.imap(func, Iter)]


def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)
        except Exception:
            error(traceback.format_exc())
            raise
        return result


class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)


def pmap(func, Iter, processes=CPUs, nice=10, *args, **kwargs):
    os.nice(nice)
    with LoggingPool(processes=processes) as P:
        results = [P.apply_async(func, (it,) + args, kwargs) for it in Iter]
        maxval = len(results)
        bar = pb.ProgressBar(
            maxval=maxval,
            widgets=[pb.Percentage(), " ", pb.Bar("=", "[", "]"), " ", pb.ETA()],
        )
        while P._cache:
            bar.update(maxval - len(P._cache))
            sleep(2)
        bar.finish()
        P.close()
    return [r.get() for r in results]


def groupby_apply(gb, func):
    keys, values = zip(*gb)
    out = pmap(func, values)
    index = gb.size().index
    if not type(out[0]) in [pd.Series, pd.DataFrame, pd.Panel]:
        return pd.DataFrame(out, index=index)
    elif len(set(index.names) - set(out[0].index.names)) == 0:
        return pd.concat(out)
    return pd.concat(out, keys=keys, names=index.names)


class _mid_fastq_iter(object):
    def __init__(self, filename, seq_id, start, stop):
        self.filename = filename
        self.seq_id = seq_id
        self.start = start
        self.stop = stop

    def __iter__(self):
        self.f = self.filename.open("rb")
        self.f.seek(self.start)
        # Find the beginning of the next FASTQ header
        lines = []
        for i, line in zip(range(5), self.f):
            lines.append(line)
            if line.startswith(self.seq_id):
                break
        else:
            raise RuntimeError(
                "Took more than 4 lines to find header in middle of FASTQ file (start pos: {:}, stop pos: {:}, file length: {:}):\n".format(
                    self.start, self.stop, self.filename.stat().st_size
                )
                + "\n".join(lines)
            )
        self.f.seek(self.f.tell() - len(line))
        return self

    def __next__(self):
        if self.f.tell() < self.stop:
            header = self.f.readline()
            dna = self.f.readline()
            self.f.readline()
            qc = self.f.readline()
            return header, dna, qc
        else:
            self.f.close()
            raise StopIteration

    def __exit__(self):
        self.f.close()


compressions = dict(gzip="gunzip", gz="gunzip", bz2="bunzip2", lzma="unxz", xz="unxz")


def fastq_map_sum(
    in_fastq,
    out_filenames,
    func,
    CPUs=CPUs - 1,
    chunks=1000,
    temp_dir_prefix="tmp",
    uncompress_input=True,
):
    """Asynchronously processes an input fastq file.

    Processes reads from a single FASTQ file by distributing the analysis work *and*
    I/O. A thread is devoted to reading the input FASTQ, threads are devoted to
    writing to every output file, and a number of worker threads communicate with
    these I/O threads to process the reads. Job queues mediate thread interactions.

    Inputs:
    -------
    in_fastq : Input FASTQ filename. Must be uncompressed.

    out_files : List of output filenames. Compression must be smart_open compliant.

    func : func(fastq_read_iter, out_filenames) -> tuple of sum-able objects.
    """
    in_file = Path(in_fastq)

    if in_file.suffix in compressions and uncompress_input:
        algorithm = compressions[in_file.suffix]
        print(
            "Uncompressing",
            in_fastq,
            "with",
            algorithm,
            "(Cannot parallel-process a compressed file)...",
        )
        os.system([algorithm, in_fastq])
        in_file = Path(in_file.stem)

    file_length = in_file.stat().st_size

    start_positions = np.linspace(0, file_length, chunks, False).round().astype(int)
    stop_positions = np.r_[start_positions[1:], file_length]

    sample = in_file.name.partition(".fastq")[0]

    Dir = Path(temp_dir_prefix + sample)

    with in_file.open("rb") as f:
        seq_id = f.readline().partition(":".encode("ascii"))[0]

    Iter = [
        (
            _mid_fastq_iter(in_file, seq_id, start, stop),
            [
                str(Dir / head / (str(start) + tail))
                for head, tail in map(os.path.split, out_filenames)
            ],
        )
        for start, stop in zip(start_positions, stop_positions)
    ]

    for head, tail in map(os.path.split, out_filenames):
        os.makedirs(head, exist_ok=True)
        (Dir / head).mkdir(parents=True, exist_ok=True)

    with multiprocessing.Pool(processes=CPUs) as P:
        rs = P.starmap_async(func, Iter, chunksize=1)
        outputs = rs.get()
    intermediate_files = list(zip(*[outfiles for it, outfiles in Iter]))
    for f, i_files in zip(out_filenames, intermediate_files):
        os.system("cat " + " ".join(i_files) + " >> " + f)
        for File in i_files:
            os.remove(File)
        (Dir / f).parent.rmdir()
    Dir.rmdir()
    return list(map(sum, zip(*outputs))) if type(outputs[0]) == tuple else sum(outputs)
