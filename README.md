# Remixavier

This repository contains code for correcting timing and channel distortion across audio signals with content in common.
In other words, it can take two related audio files, align them in time, and (approximately) correct for any difference in channel.
The python script ``remixavier.py`` is a GUI app which faciliates the process of extracing sources from a mixed audio signal given the sources which should be removed.
This setting often arises when a musician releases only an instrumental or a cappella mix of a song, but not both, and you want to either remove or isolate the vocals using the provided a cappella or instrumental mix (respectively).

### Reference

Full details of the algorithm implemented in this code are available in

C. Raffel and D. P. W. Ellis, ["Estimating Timing and Channel Distortion Across Related Signals"](http://colinraffel.com/publications/icassp2014estimating.pdf), Proceedings of the 2014 IEEE International Conference on Acoustics, Speech and Signal Processing, 2014.

If you use this code in any academic work, please cite the above paper.

The python script ``experiments.py`` contains all of the code required to generate the figures in the paper.  However, the data is not included in this github repository due to its size (>1 GB).  If you are interested in the data used in the paper, please get in touch with the authors.

### Dependencies

* [Scipy/Numpy](http://www.scipy.org/)
* [librosa](https://github.com/bmcfee/librosa)
* [mir_eval](https://github.com/craffel/mir_eval)
* [PyQt4](http://www.riverbankcomputing.co.uk/software/pyqt/download) (only required by ``remixavier.py``)
* [matplotlib](http://matplotlib.org/) (only required by ``experiments.py``)
