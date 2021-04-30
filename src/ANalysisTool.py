import sys
import os
import numpy as np
import datetime as dt

class toolBox:
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    # As this is a simple collection of functions no constructor should be 
    # needed.
    def __init__(self):
        pass

    def get_intervalSlice(dates, sliceInterval='month'):
    # def get_intervalSlice(dates, dumpInterval, sliceInterval='month'):
        ''' This functions calculates interval slices of a given time-series

        This function calculates interval slices of a given time-series, 
        aiming to mimic pythons default slice notation `array[start:stop:step]`
        to also enable multi dimensional slicing. 
        The calculation takes the models dumpintervall into account, wherefore 
        this function is working for (nearly) every time-resolution.
        The calculation is based on datetime-objects and therefore really 
        calculates the slice of individual interval, even if the time-series 
        does not start or end at the first or last of a given interval.

        Use case:
        You got a time-series of hourly data points and want to calculate the
        monthly mean. This functions returns the slices needed to divide your
        data into monthly peaces.

        Return value:
        -------------
        Slices: list
            A list with length of found intervals with related slices

        This function assumes:
        ----------------------
        -) The passed time-series covers at least on intervall!
        -) At least hourly steps / dumpIntervals! 
        -) No seconds in dates - model output is at least on full minutes!

        Parameters
        ----------
        dates : NDarray 
            the time-series as datetime-object. 
        sliceInterval: str
            defining the interval. Supported are: 'day', 'month'
        '''

        # Check is passed sliceInterval is supported and exit if not.
        supportedSliceInterval = ['day', 'month']
        if sliceInterval not in supportedSliceInterval:
            print(f'ERROR: sliceInterval "{sliceInterval}" is not supported.')
            print(f'---    supported values are: {supportedSliceInterval}')
            print(f'---    EXIT program')
            return False

        print(f'########################################################################')
        print(f'#### calculate dumpInterval and check if equal for all data-points')
        print(f'########################################################################')
        # Calculating dumpInterval
        tmp_dumpInterval = np.diff(dates, n=1)
        # Check if dumpInterval equal for all data-points
        # and if so define set dumpInterval 
        if not np.all(tmp_dumpInterval == tmp_dumpInterval[0]):
            print('ERROR: calculated dumpInterval is not equal for all data-points')
            # In case of error: break function
            return False
        else:
            dumpInterval = tmp_dumpInterval[0]
            print(f'DONE')        
        
        # Finding the first time-step of a interval, by looping over the time-series and 
        # check for each time-step if this is the first of the interval (break loop if found).
        print(f'########################################################################')
        print(f'#### Finding first of month')
        print(f'########################################################################')
        Offset = 0
        while True:
            # verify first date is first of interval.
            # idea explained with sliceInterval='month': 
            # first step - dumpInterval (or 0.5*dumpInterval) is first of month 00UTC
            # By checking 1*dumpInterval and 0.5*dumpInterval we do take into account, that
            # the time-axis of averaged model output might got shifted in between the 
            # time-bounds.
            print(f'Offset: {Offset}')
            tmp_first = dates[Offset]
            if sliceInterval == 'day':
                # using midnight as reference
                tmp_first_ref = tmp_first.replace(hour=0, minute=0, second=0)
            elif sliceInterval == 'month':
                # using the first of current month at midnight as reference
                tmp_first_ref = tmp_first.replace(day=1, hour=0, minute=0, second=0)

            if (tmp_first - dumpInterval) == tmp_first_ref:
                print(f'check step {Offset} is first of a month at midnight')
                break
            elif (tmp_first - (0.5*dumpInterval)) == tmp_first_ref:
                print(f'check step {Offset} is first of a month at midnight')
                break
            else:
                print(f'ERROR: step {Offset} is not first step of month at midnight!')
                print(f'step: {tmp_first}')
                Offset += 1
                # hard break of loop if no 'beginning' is found after 50 time-steps
                if Offset > 50:
                    break
                continue

        # Calculating the slices by looping over all time-steps and check if the 
        # next time-step belong to the next interval.
        # NWR 20210422:
        # there should be some clever and short solution for below loop ...!
        # now that we know the dumpInterval and first step is first of month...
        print(f'########################################################################')
        print(f'#### getting month series / slices')
        print(f'########################################################################')
        t_lower = Offset
        Slices = []
        for t in range(Offset, dates.shape[0]):
            # I decided to go for the solution checking current month and next month
            # to catch the case if the dateset contains one month only!
            if sliceInterval == 'day':
                currInterval = dates[t].day
                nextInterval = (dates[t] + dumpIntervall).day
            elif sliceInterval == 'month':
                currInterval = dates[t].month
                nextInterval = (dates[t] + dumpIntervall).month
            if nextInterval != currInterval:
                Slices.append(slice(t_lower,t+1,None))
                t_lower = t+1

        return Slices 