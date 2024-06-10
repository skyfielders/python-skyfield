# coding: utf-8
# jpl_mdasegment.py
"""A supporting module for for skyfield to handle SPK data type 1 (Modified
Difference Arrays) and type 21 (Extended Modified Difference Arrays)

You can get SPK files for many solar system small bodies from HORIZONS
system of NASA/JPL. See https://ssd.jpl.nasa.gov/horizons/

At the point of June. 2024, HORIZONS system creates files of type 21 for
binary SPK files by default.
You can get type 1 binary SPK file for celestial small bodies through TELNET
interface by answering back '1' for 'SPK file format'.

Author: Tontyna

This module has been developed based on Shushi Uetsuki's (whiskie14142) translation
of the FORTRAN source of the SPICE Toolkit of NASA/JPL/NAIF.
Its a extract and a combination of his Python modules spktype21 and spktype01.

skyfield : https://pypi.org/project/skyfield/
jplephem : https://pypi.org/project/jplephem/
spktype21 : https://pypi.org/project/spktype21/
spktype01 : https://pypi.org/project/spktype01/
SPICE Toolkit for FORTRAN : http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
SPK Required Reading : http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html

"""

from numpy import array, zeros, reshape, searchsorted

from jplephem.descriptorlib import reify
from jplephem.spk import BaseSegment, T0, S_PER_DAY
from jplephem.exceptions import OutOfRangeError

class MDASegment(BaseSegment):
    """Class for SPK kernel to handle data types 1 and 21 (Modified Difference Arrays)
    """
    def __init__(self, daf, source, descriptor):
        super(MDASegment, self).__init__(daf, source, descriptor)

        if self.data_type == 1:
            self.MAXDIM = 15
            # self.DLSIZE = 71
        elif self.data_type == 21:
            # type 21 has a difference line size (DLSIZE) indicating the length
            # of the MDA records in the range of
            #     [71 : (4*MAXTRM) + 1]

            # Included from 'spk21.inc' on the FORTRAN source 'spke21.f'
            MAXTRM = 25

            # 'SPK Required Reading' indicates that the penultimate element of the segment
            # is the difference line size (DLSIZE), but actual data contains there a MAXDIM.
            self.MAXDIM = int(self.daf.map_array(self.end_i - 1, self.end_i - 1))
            # self.DLSIZE = 4 * self.MAXDIM + 11

            # Q: What's the reason for limiting the table dimension to 25?
            # A: FORTRAN? History? Memory usage?
            if self.MAXDIM > MAXTRM:
                mes = ('SPKE21 \nThe input record has a maximum table dimension ' +
                    'of {0}, while the maximum supported by this routine is {1}. ' +
                    'It is possible that this problem is due to your software ' +
                    'beeing out of date.').format(self.MAXDIM, MAXTRM)
                raise RuntimeError(mes)
        else:
            raise ValueError('this class only supports SPK data types 1 and 21, ' +
                            'actual type is {0}'.format(self.data_type))

        # 4 volle bahnen plus 11 extras für tdobles (?) für die diversen
        self.DLSIZE = 4 * self.MAXDIM + 11

        # TODO: Wie kann man rauskriegen, wie weit das höchste upper_boundary reicht?
        """
        Dumm dumm dumm: self.end_second ist nicht die Zeit ab der dann nix mehr gejt!

        self.end_second ist in der Tat die letzte epoche für die noch was gerechnet werden kann
        aka epoch[recordcount-1]
        wieviele posten sie uns noch liefert, also für wieviele Seundn sie daten in sich hat???

        """


    def compute(self, tdb, tdb2=0.0):
        """Compute the component values for the time `tdb` plus `tdb2`."""
        # HA! `generate` -- für wenn die tdb ein array() ist,
        # sollte für jedes tdb ne position zurückgeliefert werden...
        for position in self.generate(tdb, tdb2):
            return position

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Compute components and differentials for time `tdb` plus `tdb2`."""
        return tuple(self.generate(tdb, tdb2))

    def generate(self, tdb, tdb2):
        """Generate components and differentials for time `tdb` plus `tdb2`.
        """

        scalar = not getattr(tdb, 'shape', 0) and not getattr(tdb2, 'shape', 0)
        if scalar:
            tdb = array((tdb,))

        eval_sec = (tdb - T0)
        eval_sec = (eval_sec + tdb2) * S_PER_DAY

        if (eval_sec < self.start_second).any() or (eval_sec >= self.end_second).any():
            raise OutOfRangeError(
                'segment only covers dates %.1f through %.1f'
                % (self.start_jd, self.end_jd),
                out_of_range_times= (eval_sec < self.start_second) | (eval_sec >= self.end_second),
            )

        # get index of the first final epoch > eval_sec,
        epoch_table = self._data
        # Hint:
        # index will never be less than 0, Even with a vereveryvery small time.
        # To exclude such out-of-rang-times, we have to check for time being smaller than self.start_second
        # Since last *final* epoch in the epoch_table == self.end_second
        # no need to check for index >= len(epoch_table)
        record_index = searchsorted(epoch_table, eval_sec, side='right')

        # collect
        positions = []
        rates = []
        for sec, idx in zip(eval_sec, record_index):
            position, rate = self._evaluate_record(sec, idx)
            positions.append(position)
            rates.append(rate)

        if scalar:
            yield positions[0]
            yield rates[0]
        else:
            # convert to nparray of  x- ,y- ,z-rows
            yield array([ row for row in zip(*positions) ])
            yield array([ row for row in zip(*rates) ])


    @reify
    def _data(self):
        if self.data_type == 1:
            offset = 0
        elif self.data_type == 21:
            # accounting for the DLSIZE
            offset = 1
        else:
            raise ValueError('this class only supports SPK data types 1 and 21')

        # number of records in this segment
        entry_count = int(self.daf.read_array(self.end_i, self.end_i))
        if entry_count < 1:
            raise RuntimeError('SPK file corrupt? Segment contains no records')

        # epoch_table contains -- distinct & increasing -- *final* epochs for all records in this segment
        # since we're using numpy.searchsorted - no need for the epoch directory
        epoch_table = self.daf.map_array(self.start_i + (entry_count * self.DLSIZE),
                                         self.start_i + (entry_count * self.DLSIZE) + entry_count - 1)
        if len(epoch_table) != entry_count:
            raise RuntimeError('SPK file corrupt? Inconsitent epoch table in segment')

        return epoch_table

    def _evaluate_record(self, epoch_time, record_index):
        """Compute position and velocity from a Modified Difference Array record

        Inputs:
            epoch_time  : Epoch time to evaluate position and velocity (seconds
                          past J2000 TDB.

            record_index: index of the MDA record to evaluate

        Returns:
            position, velocity: two  numpy arrays which contains position (km)
                                and velocity (km per day)
        """

        STATE_POS = zeros(3)
        STATE_VEL = zeros(3)

        MAXDIM = self.MAXDIM

        # TODO: reify the buffers?
        FC_BUF = zeros(MAXDIM)
        FC_BUF[0] = 1.0            # ?? unused?
        WC_BUF = zeros(MAXDIM - 1)
        W_BUF = zeros(MAXDIM + 2)

        # grab the MDA record data
        RECORD = self.daf.map_array(self.start_i + ( record_index    * self.DLSIZE),
                                    self.start_i + ((record_index+1) * self.DLSIZE) - 1)

        #
        #     Unpack the contents of the MDA array.
        #
        #        Name    Dimension  Description
        #        ------  ---------  -------------------------------
        #        TL              1  Final epoch of record
        #        G_VECTOR   MAXDIM  Stepsize function vector
        #        REFPOS          3  Reference position vector
        #        REFVEL          3  Reference velocity vector
        #        DT     MAXDIM,NTE  Modified divided difference arrays
        #        KQMAX1          1  Maximum integration order plus 1
        #        KQ            NTE  Integration order array
        #
        #     For our purposes, NTE is always 3.
        #
        TL = RECORD[0]
        G_VECTOR = RECORD[1:MAXDIM + 1]

        REFPOS = zeros(3)
        REFVEL = zeros(3)
        REFPOS[0] = RECORD[MAXDIM + 1]
        REFVEL[0] = RECORD[MAXDIM + 2]
        REFPOS[1] = RECORD[MAXDIM + 3]
        REFVEL[1] = RECORD[MAXDIM + 4]
        REFPOS[2] = RECORD[MAXDIM + 5]
        REFVEL[2] = RECORD[MAXDIM + 6]

        DT = reshape(RECORD[MAXDIM + 7 : MAXDIM * 4 + 7], (MAXDIM, 3), order='F')

        KQMAX1 = int(RECORD[4 * MAXDIM +  7])
        KQ = array([0, 0, 0])
        KQ[0] =  int(RECORD[4 * MAXDIM +  8])
        KQ[1] =  int(RECORD[4 * MAXDIM +  9])
        KQ[2] =  int(RECORD[4 * MAXDIM + 10])

        DELTA = epoch_time - TL
        TP = DELTA

#
#     This is clearly collecting some kind of coefficients.
#     The problem is that we have no idea what they are...
#
#     The G_VECTOR coefficients are supposed to be some kind of step size
#     vector.
#
#     TP starts out as the delta t between the request time
#     and the time for which we last had a state in the MDL file.
#     We then change it from DELTA  by the components of the stepsize
#     vector G_VECTOR.
#
        # what about KQMAX1 being bigger than ... MAXDIM?
        for j in range(1, KQMAX1 - 1):
            # no division by zero!
            if G_VECTOR[j-1] == 0.0:
                mes = ('SPKE01\nA value of zero was found at index {0} ' +
                'of the step size vector.').format(j)
                raise RuntimeError(mes)

            FC_BUF[j] = TP / G_VECTOR[j-1]
            WC_BUF[j-1] = DELTA / G_VECTOR[j-1]
            TP = DELTA + G_VECTOR[j-1]

#
#     Collect KQMAX1 reciprocals.
#
        for j in range(1, KQMAX1 + 1):
            W_BUF[j-1] = 1.0 / float(j)
#
#     Compute the W_BUF(K) terms needed for the position interpolation
#     (Note,  it is assumed throughout this routine that KS, which
#     starts out as KQMAX1-1 (the ``maximum integration'')
#     is at least 2.
#
        JX  = 0
        KS  = KQMAX1 - 1
        KS1 = KS - 1
        while KS >= 2:
            JX += 1
            for j in range(1, JX + 1):
                W_BUF[j+KS-1] = FC_BUF[j] * W_BUF[j+KS1-1] - WC_BUF[j-1] * W_BUF[j+KS-1]
            KS = KS1
            KS1 -= 1
#
#     Perform position interpolation: (Note that KS = 1 right now.
#     We don't know much more than that.)
#
        for i in range(1, 3 + 1):
            kqq = KQ[i-1]
            sum = 0.0
            for j in range(kqq, 0, -1):
                sum += DT[j-1, i-1] * W_BUF[j+KS-1]

            STATE_POS[i-1] = REFPOS[i-1] + DELTA * (REFVEL[i-1] + DELTA * sum)
#
#     Again we need to compute the W_BUF(K) coefficients that are
#     going to be used in the velocity interpolation.
#     (Note, at this point, KS = 1, KS1 = 0.)
#
        for j in range(1, JX + 1):
            W_BUF[j+KS-1] = FC_BUF[j] * W_BUF[j+KS1-1] - WC_BUF[j-1] * W_BUF[j+KS-1]

        KS = KS - 1
#
#     Perform velocity interpolation:
#
        for i in range(1, 3 + 1):
            kqq = KQ[i-1]
            sum = 0.0

            for j in range(kqq, 0, -1):
                sum += DT[j-1, i-1] * W_BUF[j + KS-1]

            STATE_VEL[i-1] = REFVEL[i-1] + DELTA * sum
#
#     That's all folks.  We don't know why we did anything, but
#     at least we can tell structurally what we did.
#

        # MDA array gives velocities in kilometers per second,
        # jplephem returns (and skyfield expects) velocities in kilometers per day.
        return STATE_POS, STATE_VEL * S_PER_DAY
