"""
    :platform: Unix, Windows
    :synopsis: This module contains utility functions for Bertini_real, including functions:
      * to generate BRData*.pkl 
      * to enhance lists
"""

def next_filenumber():
    """ Keep track on the next filenumber for BRData

        :rtype: An integer which is the next file number of BRData*.pkl

    """
    import fnmatch
    import os

    pattern = 'BRdata*.pkl'

    files = os.listdir('.')
    highest_number = -1

    for name in files:
        if fnmatch.fnmatch(name, pattern):
            try:
                current_number = int(name[6:-4])
                if current_number > highest_number:
                    highest_number = current_number
            except ValueError:
                continue

    return highest_filenumber()+1


def highest_filenumber():
    """ Get the highest/most recent file number of BRData

        :rtype: The highest integer of BRData*.pkl

    """
    import fnmatch
    import os

    pattern = 'BRdata*.pkl'

    files = os.listdir('.')
    highest_number = -1

    for name in files:
        if fnmatch.fnmatch(name, pattern):
            try:
                current_number = int(name[6:-4])
                if current_number > highest_number:
                    highest_number = current_number
            except ValueError:
                continue

    return highest_number



class ReversableList(list):
    """ Create a ReversableList object for reversing order of data 

        :param list: The list to be read.

    """

    def reverse(self):
        return list(reversed(self))
