
import logomaker
import pandas
import numpy
from numpy.lib.stride_tricks import as_strided

from matplotlib import pyplot

colors = {
    'A': 'blue',
    'C': 'yellow',
    'G': 'green',
    'T': 'red',
    '-': 'black'
}


def compute_end_pos(read):
    return read.start + read.length


def window(
    array:numpy.ndarray, window:int,
    shift:int=1, copy:bool=False
):
    shape = (array.size - window + 1, window)
    stride = array.strides * 2

    view = as_strided(
        array, strides=stride, shape=shape
    )[0::shift]

    if copy:
        return view.copy()

    else:
        return view


def create_seqlogo_dataframe(sequences):
    seqlen = len(sequences[0])
    idx = pandas.Index(data=range(seqlen), name='pos')
    basemat = numpy.array([list(x) for x in sequences]).T
    df = pandas.DataFrame(
        numpy.zeros(shape=(seqlen, len(colors))), index=idx, columns=colors.keys()
    )

    for x in colors:
        df[x] = (basemat == x).sum(axis=1)

    freqs = df / len(sequences)
    log_freqs = numpy.nan_to_num(numpy.log2(freqs))

    correction = (1 / numpy.log(2)) * (3/(2*len(sequences)))
    row_scale = numpy.log2(5) - (- (freqs * log_freqs) + correction).sum(axis=1)

    return freqs.multiply(row_scale, axis=0)


def create_logo(sequences, figsize=(10, 2.5), save=None):

    fig, ax = pyplot.subplots(1, figsize=figsize)

    df = create_seqlogo_dataframe(sequences)

    logo = logomaker.Logo(df, color_scheme=colors, ax=ax)
    logo.style_xticks(anchor=0, spacing=5, rotation=45)
    logo.ax.set_xlim([-1, len(df)])
    logo.ax.set_ylabel('information (bits)')

    return logo, fig
